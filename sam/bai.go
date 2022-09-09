package sam

import (
	"bytes"
	"encoding/binary"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"math"
	"sort"
)

// bai files are an index format for bam files
// Specs for bai can be found in section 5.2 on page 19 of:
// https://raw.githubusercontent.com/samtools/hts-specs/master/SAMv1.pdf

// magicBai is a 4 byte sequence at the start of a bai file
const magicBai string = "BAI\u0001"

// alias for ease of use and readability
var le = binary.LittleEndian

// Bai is an index of the BAM format that allows for efficient
// random access of reads given a query genomic range. The Bai
// struct should not be accessed directly, but used through
// query functions present within the sam package of gonomics.
type Bai struct {
	refs            []reference
	hasNoCoordReads bool   // true if optional field below exists
	noCoordReads    uint64 // number of unplaced unmapped reads
}

// reference sequence being binned.
type reference struct {
	sliceIdx    uint32 // index in Bai.refs
	head        bin    // a dummy bin at the top of tree to organize otherwise unlinked bins.
	bins        []bin
	binIdIdx    map[int]int // binIdIdx[bin.Id] == bin.sliceIdx
	intervalOff []uint64
	noCoord     noCoordBin
}

// bin containing all reads that are contained within a variable genomic size.
// bins are organized into a tree structure for efficient queries by genomic coordinates.
type bin struct {
	sliceIdx uint32 // index in reference.bins
	id       uint32 // distinct bin id (only unique within reference). see const for id -> size.
	refStart uint32 // genomic start coordinate of bin (0-base, closed)
	refEnd   uint32 // genomic end coordinate of bin (0-base, open)
	chunks   []chunk
	parent   *bin
	children []*bin
}

// chunk represents the start and end virtual offset of a contiguous
// block of bam records that fall within a given bin. All records within
// a bin are not necessarily contiguous, therefore a bin may have multiple chunks.
type chunk struct {
	start uint64 // virtual offset to start of chunk
	end   uint64 // virtual offset to end of chunk
}

// noCoordBin is a special 'dummy' bin with id == 37450 that stores
// metadata about the number of mapped/unmapped reads for a given reference
// as well as the start and end virtual offsets of the unmapped reads in
// the bam file.
type noCoordBin struct {
	start       uint64 // virtual offset to start of chunk
	end         uint64 // virtual offset to end of chunk
	numMapped   uint64
	numUnmapped uint64
}

// ReadBai reads an entire bai file into a bytes buffer and parses
// the buffer into a Bai struct which can be used for efficient
// random access of reads in a bam file.
func ReadBai(filename string) Bai {
	file := fileio.EasyOpen(filename)
	r := new(bytes.Buffer)
	_, err := r.ReadFrom(file) // read entire file to bytes buffer
	exception.PanicOnErr(err)
	err = file.Close() // close read file
	exception.PanicOnErr(err)

	if string(r.Next(4)) != magicBai {
		log.Fatalf("malformed bai header in '%s'", filename)
	}

	bai := parseBai(r)

	if r.Len() == 8 { // check for presence of optional field
		bai.hasNoCoordReads = true
		bai.noCoordReads = le.Uint64(r.Next(8))
	}

	if r.Len() != 0 {
		log.Fatalf("ERROR: extra '%d' bytes found in bai file.\nfile may be malformed: '%s'\n", r.Len(), filename)
	}

	for i := range bai.refs {
		bai.refs[i] = annotateBinRanges(bai.refs[i])
		bai.refs[i] = assembleTree(bai.refs[i])
		bai.refs[i].binIdIdx = make(map[int]int)
		for j := range bai.refs[i].bins {
			bai.refs[i].binIdIdx[int(bai.refs[i].bins[j].id)] = j
		}
	}
	return bai
}

// parseBai processes the bytes buffer (after reading the magic bytes)
// to generate a Bai struct.
func parseBai(r *bytes.Buffer) Bai {
	var bai Bai
	numRefs := le.Uint32(r.Next(4)) // note le is an alias for binary.LittleEndian
	bai.refs = make([]reference, numRefs)

	var i uint32
	for i = 0; i < numRefs; i++ {
		bai.refs[i] = parseReference(r)
		bai.refs[i].sliceIdx = i
	}
	return bai
}

// parseReference processes the bytes buffer to generate a reference struct.
func parseReference(r *bytes.Buffer) reference {
	var ref reference
	numBins := le.Uint32(r.Next(4)) // note le is an alias for binary.LittleEndian
	ref.bins = make([]bin, numBins)

	var i uint32
	for i = 0; i < numBins; i++ {
		ref.bins[i] = parseBin(r)
		ref.bins[i].sliceIdx = i
		//TODO parents and children

		// handle if bin is an optional noCoordBin
		if ref.bins[i].id == 37450 { // the magic bin number
			ref.noCoord.start = ref.bins[i].chunks[0].start
			ref.noCoord.end = ref.bins[i].chunks[0].end
			ref.noCoord.numMapped = ref.bins[i].chunks[1].start
			ref.noCoord.numUnmapped = ref.bins[i].chunks[1].end
			ref.bins[i] = bin{}                   // remove from bin list and replace with empty bin
			ref.bins = ref.bins[:len(ref.bins)-1] // decrease size of ref.bins by 1
			i--
			numBins--
		}
	}

	numIntervals := le.Uint32(r.Next(4))
	ref.intervalOff = make([]uint64, numIntervals)
	for i = 0; i < numIntervals; i++ {
		ref.intervalOff[i] = parseIntervalOffset(r)
	}

	return ref
}

// parseBin processes the bytes buffer to generate a bin struct.
func parseBin(r *bytes.Buffer) bin {
	var b bin
	b.id = le.Uint32(r.Next(4)) // note le is an alias for binary.LittleEndian
	numChunks := le.Uint32(r.Next(4))
	b.chunks = make([]chunk, numChunks)

	var i uint32
	for i = 0; i < numChunks; i++ {
		b.chunks[i] = parseChunk(r)
	}

	return b
}

// parseChunk processes the bytes buffer to generate a chunk struct.
func parseChunk(r *bytes.Buffer) chunk {
	var c chunk
	c.start = le.Uint64(r.Next(8)) // note le is an alias for binary.LittleEndian
	c.end = le.Uint64(r.Next(8))
	return c
}

// parseIntervalOffset processes the bytes buffer to generate an interval offset.
func parseIntervalOffset(r *bytes.Buffer) uint64 {
	return le.Uint64(r.Next(8)) // note le is an alias for binary.LittleEndian
}

// to reduce chance of errors in the following function
const million uint32 = 1000000
const thousand uint32 = 1000

// annotateBinRanges adds the start and end coordinates to each bin in a bai file.
func annotateBinRanges(r reference) reference {
	for i := range r.bins {
		switch {
		case r.bins[i].id > 37448:
			log.Panicf("bin id: '%d' overflow. bai may be malformed\n", r.bins[i].id)

		case r.bins[i].id > 4680: // size is 16kbp
			r.bins[i].refStart = (r.bins[i].id - 4681) * (16 * thousand)
			r.bins[i].refEnd = r.bins[i].refStart + (16 * thousand)

		case r.bins[i].id > 584: // size is 128kbp
			r.bins[i].refStart = (r.bins[i].id - 585) * (128 * thousand)
			r.bins[i].refEnd = r.bins[i].refStart + (128 * thousand)

		case r.bins[i].id > 72: // size is 1Mbp
			r.bins[i].refStart = (r.bins[i].id - 73) * million
			r.bins[i].refEnd = r.bins[i].refStart + million

		case r.bins[i].id > 8: // size is 8Mbp
			r.bins[i].refStart = (r.bins[i].id - 9) * (8 * million)
			r.bins[i].refEnd = r.bins[i].refStart + (8 * million)

		case r.bins[i].id > 0: // size is 64Mbp
			r.bins[i].refStart = (r.bins[i].id - 1) * (64 * million)
			r.bins[i].refEnd = r.bins[i].refStart + (64 * million)

		case r.bins[i].id == 0: // size is 512Mbp
			r.bins[i].refStart = 0
			r.bins[i].refEnd = 512 * million

		}
	}
	return r
}

// assembleTree arranges bins into an R tree structure.
func assembleTree(r reference) reference {
	// establish a slice of bins that is sorted by genomic position (small -> large)
	// and secondly sorted by size (large -> small)
	proxy := make([]*bin, len(r.bins)) // sortable proxy slice
	for i := range r.bins {
		proxy[i] = &r.bins[i]
	}
	sort.Slice(proxy, func(i, j int) bool {
		switch {
		case proxy[i].refStart < proxy[j].refStart:
			return true

		case proxy[i].refStart > proxy[j].refStart:
			return false

		case proxy[i].size() > proxy[j].size():
			return true

		default:
			return false
		}
	})

	r.head = bin{refEnd: math.MaxUint32} // allocate dummy head bin

	// the parents slice will keep track of the last node encountered
	// prior to switching bin sizes. When going down in bin size a new
	// bin is appended to the slice. When going up in bin size, the last
	// parent node is trimmed.
	var parents []*bin
	var currParent *bin

	// the head node is set to have a larger size than any valid bin
	// therefore it will always be the first bin in parents and will
	// not be removed.
	var prev *bin = &r.head

	for i := range proxy {
		switch {
		// if prev.size() == b.size() no bin size change so parents unchanged

		case prev.size() > proxy[i].size(): // bin size decreased. parents grows
			parents = append(parents, prev)      // add prev to parent slice
			currParent = parents[len(parents)-1] // set curr parent

		case prev.size() < proxy[i].size(): // bin size increased. parents shrinks
			// any node without a direct parent, should point to the head node
			// so we should never trim the head node from the parents slice.
			for len(parents) > 1 {
				parents = parents[:len(parents)-1]   // trim last parent
				currParent = parents[len(parents)-1] // set curr parent

				// The following statement handles cases where bin size increases
				// by multiple increments. e.g. if we went from bin size of
				// 64Mbp -> 8Mbp -> 128kbp -> 512Mbp, the parents slice arriving at
				// the final position would be: [head, 64bin, 8bin]. Without the below
				// line (and the for loop) we would calculate that the parent of the
				// 512Mbp bin is the 64Mbp bin (as 8bin gets trimmed) which is incorrect.
				// Instead we continue to trim by 1 until the parent size > current bin size.
				if currParent.size() > proxy[i].size() {
					break
				}
			}
		}
		proxy[i].parent = currParent
		currParent.children = append(currParent.children, proxy[i])
		prev = proxy[i]
	}
	return r
}

// size returns the genomic size a bin.
func (b bin) size() uint32 {
	return b.refEnd - b.refStart
}
