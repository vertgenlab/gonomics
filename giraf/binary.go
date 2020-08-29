package giraf

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dnaThreeBit"
)

// Compressed with BGZF with TOC storing reads in blocks

type binaryGirafHeader struct {
	magic          [4]byte  // magic string to verify file is a binary giraf
	refChecksum    [16]byte // md5 checksum to make sure correct gg file is used for decompression
	textLen        uint32   // length of below text field
	text           string   // plain header text in giraf file
	qNamePrefixLen uint32   // length of below qNamePrefix field
	qNamePrefix    string   // prefix of read name, may be updated to algorithmic compression of read names
}

type binaryGiraf struct {
	blockSize     uint32 // total length of the alignment record, excluding this field
	hasNamePrefix bool   // bool to determine whether qNamePrefix is valid for this read
	qNameLen      uint8  // length of below qName field
	qName         string // name of read, qNamePrefix is prepended if hasNamePrefix
	//Question: can we drop QStart and QEnd? Do these serve a different purpose from TStart and TEnd in the Path?
	flag uint8 // bitwise flags
	// PosStrand is dropped, stored in flag
	tStart      uint32               // start position on first node (nodes[0]) in path, zero base
	tEnd        uint32               // end position on last node (nodes[len(nodes)-1]) in path, zero base
	pathLen     uint16               // length of below path
	path        []uint32             // list of nodes in alignment, len of path is = pathLen
	numCigarOps uint16               // number of cigar operations
	byteCigar   []cigar.ByteCigar    // run-length encoding of alignment. Can be copied directly from Giraf struct IF the aligner encodes mismatches into the cigar.
	fancySeq    dnaThreeBit.ThreeBit // three bit encoding of 'fancy-cigar' seq determined at time of compression. Length will be determined by Op's in byteCigar
	alnScore    int64                // alignment quality score
	mapQ        uint8                // mapping quality score
	// Seq field is dropped, can be determined from aln
	readLen uint32        // length of read
	qual    []cigar.ByteCigar       // phred-scaled base qualities. run-length encoding using ByteCigar
	notes   []binaryNotes // the notes field will be identical to the SAM notes field. See section 4.2.4 in SAM specs for details
}

type binaryNotes struct {
	tag     [2]byte // tag ID
	tagType byte    // type of data encoded. See section 4.2.4 in SAM specs
	data    []byte  // data attributed to the the tag. See section 4.2.4 in SAM specs for how to decode
}
