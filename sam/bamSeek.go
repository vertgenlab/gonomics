package sam

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"log"
	"sort"
)

// SeekBamRegion returns a slice of reads that overlap the input region. SeekBamRegion will advance the
// reader as necessary so beware continuing linear reading after a call to SeekBamRegion.
func SeekBamRegion(br *BamReader, bai Bai, chrom string, start, end uint32) []Sam {
	return SeekBamRegionRecycle(br, bai, chrom, start, end, nil)
}

// SeekBamRegionRecycle functions identially to SeekBamRegion, but inputs a slice of Sam structs that
// which will be reused. Avoids excessive memory allocations on repeat seek calls.
func SeekBamRegionRecycle(br *BamReader, bai Bai, chrom string, start, end uint32, recycledSams []Sam) []Sam {
	if start > end {
		log.Panicf("ERROR: SeekBamRegion input start > end. %d > %d\n", start, end)
	}
	var err error
	var coffset, uoffset, cEndOffset, linearIndexMinCOffset uint64
	var ans []Sam = recycledSams[0:cap(recycledSams)] // expand recycledSams to cap
	var ansIdx int = 0

	refIdx := chromInfo.SliceToMap(br.refs)[chrom].Order
	bins := regionToBins(int(start), int(end))
	ref := bai.refs[refIdx]
	var i, j, binIdx int
	var ok bool

	//For an alignment starting beyond 64Mbp, we always need to seek to some chunks in bin 0, which can be
	//avoided by using a linear index. In the linear index, for each tiling 16384bp window on the reference, we
	//record the smallest file offset of the alignments that overlap with the window. Given a region [rbeg, rend), we
	//only need to visit a chunk whose end file offset is larger than the file offset of the 16kbp window containing
	//rbeg.
	//	With both binning and linear indices, we can retrieve alignments in most of regions with just one seek
	//call.
	linearIndexMinCOffset = ref.intervalOff[start/16384] >> 16

	var curr *Sam
	var madeNewSam bool
	for i = range bins { // retrieve bins that may contain overlapping reads
		if _, ok = ref.binIdIdx[bins[i]]; !ok {
			continue
		}
		binIdx = ref.binIdIdx[bins[i]]
		for j = range ref.bins[binIdx].chunks { // check all chunks with each bin
			uoffset = ref.bins[binIdx].chunks[j].start & 0xFFFF // byte offset into uncompressed data stream
			coffset = ref.bins[binIdx].chunks[j].start >> 16    // byte offset from start to bgzf block
			//uEndOffset = ref.bins[binIdx].chunks[j].end & 0xFFFF // TODO minor efficiency benefit by including this in linear index calculation
			cEndOffset = ref.bins[binIdx].chunks[j].end >> 16
			if cEndOffset < linearIndexMinCOffset {
				continue
			}
			_, err = br.zr.Seek(int64(coffset), io.SeekStart)
			exception.PanicOnErr(err)
			err = br.zr.ReadBlock(br.blk)
			exception.PanicOnErr(err)
			br.blk.Next(int(uoffset)) // advance in block
			br.intermediate.Reset()
			for { // read sam records and append to answer until past the query region
				madeNewSam = false
				// check to see if we have any recycled structs available
				if ansIdx < len(ans) {
					curr = &ans[ansIdx]
				} else { // if none available make a new one
					curr = new(Sam)
					madeNewSam = true
				}
				_, err = DecodeBam(br, curr)
				if err == io.EOF {
					break
				}
				if curr.RName == chrom && curr.GetChromEnd() > int(start) && curr.GetChromStart() < int(end) {
					if madeNewSam { // if we used a recycled sam, it is already is the ans slice. Only need to append if it is new
						ans = append(ans, *curr)
					}
					ansIdx++
				}
				if (curr.RName == chrom && curr.GetChromStart() >= int(end)) || curr.RName != chrom {
					break
				}
			}
		}
	}
	ans = ans[:ansIdx]
	ans = deduplicate(ans) // TODO to improve efficiency this deduplication should be removed in favor of tracking the
	// current file offset in DecodeBam and breaking after passing cEndOffset and uEndOffset.
	// We are currently reading past the end of the designated chunk; this gives duplicate values
	// but gives a correct answer after deduplication. By stopping the read at uEndOffset we can
	// avoid deduplication.
	return ans
}

// deduplicate removes duplicated sam records in a slice
func deduplicate(s []Sam) []Sam {
	ans := make([]Sam, 0, len(s))
	sort.Slice(s, func(i, j int) bool {
		switch {
		case s[i].QName < s[j].QName:
			return true
		case s[i].QName > s[j].QName:
			return false
		default: // names match but could be pairs, return fwd < rev
			return IsForwardRead(s[i])
		}
	})
	for i := range s {
		if len(ans) == 0 || !(s[i].QName == ans[len(ans)-1].QName && IsForwardRead(s[i]) == IsForwardRead(ans[len(ans)-1])) {
			ans = append(ans, s[i])
		}
	}
	return ans
}

// regionToBins returns all bins in which reads overlapping the input region may exist.
// Adapted from C code in SAM file specifications.
func regionToBins(beg, end int) []int {
	if beg == -1 && end == 0 {
		return []int{0, 0, 8, 72, 584, 4680}
	}
	var ans []int
	var k int
	end--
	for k = 1 + (beg >> 26); k <= 1+(end>>26); k++ {
		ans = append(ans, k)
	}
	for k = 9 + (beg >> 23); k <= 9+(end>>23); k++ {
		ans = append(ans, k)
	}
	for k = 73 + (beg >> 20); k <= 73+(end>>20); k++ {
		ans = append(ans, k)
	}
	for k = 585 + (beg >> 17); k <= 585+(end>>17); k++ {
		ans = append(ans, k)
	}
	for k = 4681 + (beg >> 14); k <= 4681+(end>>14); k++ {
		ans = append(ans, k)
	}
	return ans
}
