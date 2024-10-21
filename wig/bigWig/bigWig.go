package bigWig

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
	DataCount         uint32            // The number of sections in the bigWig.
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

// BinaryWigSectionHeader contains information stored in the header of each Wig segment in a bigWig file.
type BinaryWigSectionHeader struct {
	ChromId    uint32 // Numerical ID for chromosome/contig
	ChromStart uint32 // Start position (0 based)
	ChromEnd   uint32 // end of item. Equal to ChromStart + ItemSize
	ItemStep   uint32 // Spaces between start of adjacent items in fixedStep sections.
	ItemSpan   uint32 // Number of bases in item in fixedStep and varStep sections.
	Type       uint8  // Section type. 1 for bedGraph, 2 for varStep, 3 for fixedStep.
	Reserved   uint8  // Currently 0.
	ItemCount  uint16 // Number of items in section.
}
