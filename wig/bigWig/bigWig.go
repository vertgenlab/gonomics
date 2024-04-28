package bigWig

import (
	"encoding/binary"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

// Note that package bigWig is written to parse bigWig format version 4.
// specs for version 4 are found in this file: https://github.com/ucscGenomeBrowser/kent/blob/master/src/inc/bbiFile.h
// the official documentation for the bigWig file format (version 3)
// can be found in the supplementary materials of the following
// publication:
// Kent, W.J., et al. BigWig and BigBed: enabling browsing of large distributed
//datasets. (2010) Bioinformatics, 26(17) 2204-2207.

// these are the magic number for bigWig files
const bigWigMagic = 2291137574             // little endian
const bigWigMagicBigEndian = 654086024     // big endian, not supported by this package
const chromTreeMagic = 2026540177          // little endian
const chromTreeMagicBigEndian = 2441923192 // big endian is not supported by this package

// BigWig represents the data of a BigWig file in a data structure.
type BigWig struct {
	BbiHeader         BbiHeader         // Contains high-level information about file and offsets to various parts of the file.
	ZoomHeaders       []ZoomHeader      // One for each zoom level built into the file, as specified in the BbiHeader.
	TotalSummaryBlock TotalSummaryBlock // Statistical summary of the entire file.
	ChromTreeHeader   ChromTreeHeader   // Index of chromosomes, their sizes, unique IDs. B+ tree data structure.
	ChromTreeNodes    []ChromTreeNode   // Nodes of the chromosome tree.
}

// BbiHeader represents the bbi header of the bigWig data file and contains metadata and offset values for file navigation.
type BbiHeader struct {
	Magic                uint32 // magic bytes at the beginning of file. Must match constant above to be a valid bigWig.
	Version              uint16 // Specifies the wig version. This package was built following specs for bigWig version 4.
	ZoomLevels           uint16 // Number of zoom levels built into the file.
	ChromosomeTreeOffset uint64 // Offset in file to chromosome B+ tree index.
	FullDataOffset       uint64 // Offset to the main dta. Points to the dataCount.
	FullIndexOffset      uint64 // Offset to R tree index of items.
	FieldCount           uint16 // Number of fields in a BED file. For bigWig, this should always be zero.
	DefinedFieldCount    uint16 // Number of fields that are predefined BED files. this should also be zero for bigWig files.
	AutoSqlOffset        uint64 // From specs: Offset to zero-terminated string with .as spec. Used for bigBed, not used in bigWig.
	TotalSummaryOffset   uint64 // Offset to overall file summary data block.
	UncompressBufferSize uint32 // Maximum size of decompression buffer needed (nonzero on compressed files).
	ExtensionOffset      uint64 // Offset to header extension 0 if no such extension.
}

// ZoomHeader immediately follows the BbiHeader in the file, one for each of the ZoomLevels, as specified in the BbiHeader.
type ZoomHeader struct {
	ReductionLevel uint32 // Number of bases summarized in each reduction level.
	Reserved       uint32 // Reserved for future expansion. 0 in bigWig.
	DataOffset     uint64 // Position of zoomed data in file.
	IndexOffset    uint64 // Position of zoomed data index in file.
}

// TotalSummaryBlock provides an overall statistical summary of the file contents.
// Mean and standard deviation values can be quickly calculated from these values.
type TotalSummaryBlock struct {
	BasesCovered uint64  // Number of bases for which there is data.
	MinVal       float64 // Minimum value in the file.
	MaxVal       float64 // Maximum value in the file.
	SumData      float64 // Sum of all values in the file.
	SumSquares   float64 // Sum of all squares of values in the file.
}

// ChromTreeHeader specifies the header for the Chromosome B+ tree header. It starts at the offset specified by
// chromosomeTreeOffset in the common header.
// This struct contains information about the structure of the chromosome tree nodes, which immediately follow this header
// in the file.
type ChromTreeHeader struct {
	Magic     uint32 // expects 2026540177 (Little endian). Package throws error if big endian header is found.
	BlockSize uint32 // Number of children per block (not byte size of the block).
	KeySize   uint32 // Number of significant bytes in key. Minimum prefix size needed to distinguish one chromosome name from another.
	ValSize   uint32 // Size of the value being indexed. Set to 8.
	ItemCount uint64 // The number of chromosomes / contigs.
	Reserved  uint64 // Reserved for future expansions. Must be 0 in bigWig.
}

// ChromTreeNode specifies the structure of nodes of the chromosome B+ tree. The first node immediately
// follows the ChromTreeHeader in the binary file.
type ChromTreeNode struct {
	IsLeaf   bool            // true if leaf, false if non-leaf.
	Reserved byte            // Reserved for future format expansion. Must be 0 in valid bigWig files.
	Count    uint16          // Number of ChromTreeItems in the node.
	Items    []ChromTreeItem // This contains all the items in the current node
}

// ChromTreeItem contains information stored within a ChromTreeNode. Note that nodes are
// either Leaf or Non-Leaf nodes, and some fields are used in only one type of node.
type ChromTreeItem struct {
	Key         []byte // This stores the first few characters of the chromosome name, padded with zeros where necessary. The length of this slice is defined by the KeySize value in ChromTreeHeader.
	ChromId     uint32 // Leaf Only. Numerical ID for chromosome / contig.
	ChromSize   uint32 // Leaf Only.Number of bases in chromosome / contig.
	ChildOffset uint64 // Non-Leaf Only. Offset to child node.
}

// Read parses a BigWig struct from an input filename.
func Read(filename string) BigWig {
	var err error
	var answer = BigWig{}
	file := fileio.EasyOpen(filename)
	answer.BbiHeader = readBbiHeader(file)
	for currZoomLevel := 0; currZoomLevel < int(answer.BbiHeader.ZoomLevels); currZoomLevel++ {
		answer.ZoomHeaders = append(answer.ZoomHeaders, readZoomHeader(file))
	}
	answer.TotalSummaryBlock = readTotalSummary(file)
	answer.ChromTreeHeader = readChromTreeHeader(file)
	answer.ChromTreeNodes = readChromTreeNodes(file, answer.ChromTreeHeader)
	err = file.Close()
	exception.PanicOnErr(err)
	return answer
}

// We currently support only reading bbi headers in LittleEndian.
func readBbiHeader(file *fileio.EasyReader) BbiHeader {
	var header = BbiHeader{}
	err := binary.Read(file, binary.LittleEndian, &header.Magic)
	exception.PanicOnErr(err)
	if header.Magic == bigWigMagicBigEndian {
		log.Fatalf("Error: bigWig file appears to be in big endian. Current functionality only supports little endian bigWig files.\n")
	}
	if header.Magic != bigWigMagic {
		log.Fatalf("Error: bigWig magic was not as expected. Found: %v. Expected: %v.\n", header.Magic, bigWigMagic)
	}
	err = binary.Read(file, binary.LittleEndian, &header.Version)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.ZoomLevels)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.ChromosomeTreeOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.FullDataOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.FullIndexOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.FieldCount)
	exception.PanicOnErr(err)
	if header.FieldCount != 0 {
		log.Fatalf("Error: bigWig header fieldCount field must be zero. Found: %v.\n", header.FieldCount)
	}
	err = binary.Read(file, binary.LittleEndian, &header.DefinedFieldCount)
	exception.PanicOnErr(err)
	if header.DefinedFieldCount != 0 {
		log.Fatalf("Error: bigWig header definedFieldCount field must be zero. Found: %v.\n", header.DefinedFieldCount)
	}
	err = binary.Read(file, binary.LittleEndian, &header.AutoSqlOffset)
	exception.PanicOnErr(err)
	if header.AutoSqlOffset != 0 {
		log.Fatalf("Error: bigWig header AutoSeqlOffset field must be zero. Found: %v.\n", header.AutoSqlOffset)
	}
	err = binary.Read(file, binary.LittleEndian, &header.TotalSummaryOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.UncompressBufferSize)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.ExtensionOffset)
	exception.PanicOnErr(err)
	return header
}

func readZoomHeader(file *fileio.EasyReader) ZoomHeader {
	var answer = ZoomHeader{}
	var err error
	err = binary.Read(file, binary.LittleEndian, &answer.ReductionLevel)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.Reserved)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.DataOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.IndexOffset)
	exception.PanicOnErr(err)
	return answer
}

func readTotalSummary(file *fileio.EasyReader) TotalSummaryBlock {
	var answer = TotalSummaryBlock{}
	var err error
	err = binary.Read(file, binary.LittleEndian, &answer.BasesCovered)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.MinVal)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.MaxVal)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.SumData)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.SumSquares)
	exception.PanicOnErr(err)
	return answer
}

func readChromTreeHeader(file *fileio.EasyReader) ChromTreeHeader {
	var answer = ChromTreeHeader{}
	var err error
	err = binary.Read(file, binary.LittleEndian, &answer.Magic)
	exception.PanicOnErr(err)
	if answer.Magic == chromTreeMagicBigEndian {
		log.Fatalf("Error: bigWig package found a big endian chromosome tree header. Currently our package only supports little endian bigWig files.\n")
	}
	if answer.Magic != chromTreeMagic {
		log.Fatalf("Error: Expected to find chromosome tree magic number (2026540177) in file. Found: %v.\n", answer.Magic)
	}
	err = binary.Read(file, binary.LittleEndian, &answer.BlockSize)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.KeySize)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.ValSize)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.ItemCount)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.Reserved)
	exception.PanicOnErr(err)
	if answer.Reserved != 0 {
		log.Fatalf("Error: expected chromosome tree header reserved field to be equal to 0. Found: %v.\n", answer.Reserved)
	}
	return answer
}

func readChromTreeNodes(file *fileio.EasyReader, header ChromTreeHeader) []ChromTreeNode {
	var answer = make([]ChromTreeNode, 0)
	var err error
	var countIdx uint16
	var currByte uint32
	var nodeIdx uint64
	var currNode ChromTreeNode
	var currItem ChromTreeItem
	for nodeIdx = 0; nodeIdx < header.ItemCount; nodeIdx++ {
		currNode = ChromTreeNode{}
		err = binary.Read(file, binary.LittleEndian, &currNode.IsLeaf)
		exception.PanicOnErr(err)
		err = binary.Read(file, binary.LittleEndian, &currNode.Reserved)
		exception.PanicOnErr(err)
		err = binary.Read(file, binary.LittleEndian, &currNode.Count)
		exception.PanicOnErr(err)

		for countIdx = 0; countIdx < currNode.Count; countIdx++ {
			currItem = ChromTreeItem{}
			// first we parse the key (same for leaf and non-leaf nodes, and comes first)
			currItem.Key = make([]byte, header.KeySize)
			for currByte = 0; currByte < header.KeySize; currByte++ {
				err = binary.Read(file, binary.LittleEndian, &currItem.Key[currByte])
				exception.PanicOnErr(err)
			}
			if currNode.IsLeaf {
				err = binary.Read(file, binary.LittleEndian, &currItem.ChromId)
				exception.PanicOnErr(err)
				err = binary.Read(file, binary.LittleEndian, &currItem.ChromSize)
				exception.PanicOnErr(err)
			} else {
				err = binary.Read(file, binary.LittleEndian, &currItem.ChildOffset)
				exception.PanicOnErr(err)
			}
			currNode.Items = append(currNode.Items, currItem)
		}

		answer = append(answer, currNode)
	}
	return answer
}
