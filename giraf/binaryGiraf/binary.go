package binaryGiraf

// Package binaryGiraf provides functions for reading and writing binary giraf files

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dnaThreeBit"
)

// The structs in this file are only present for documentation and should not be used

// Compressed with BGZF with TOC storing reads in blocks

type binGirafHeader struct {
	magic          [4]byte  // magic string to verify file is a binary giraf
	refChecksum    [16]byte // md5 checksum to make sure correct gg file is used for decompression
	textLen        uint32   // length of below text field
	text           string   // plain header text in giraf file
	qNamePrefixLen uint32   // length of below qNamePrefix field
	qNamePrefix    string   // prefix of read name, may be updated to algorithmic compression of read names
}

type binGiraf struct {
	blockSize uint32 // total length of the alignment record, excluding this field
	// THIS HAS BEEN MOVED TO FLAG //hasNamePrefix bool   // bool to determine whether qNamePrefix is valid for this read
	qNameLen uint8  // length of below qName field
	qName    string // name of read, qNamePrefix is prepended if hasNamePrefix
	// QStart and QEnd are determined by clipped bases in Cigar
	flag uint8 // bitwise flags
	// PosStrand is dropped, stored in flag
	tStart      uint32               // start position on first node (nodes[0]) in path, zero base
	tEnd        uint32               // end position on last node (nodes[len(nodes)-1]) in path, zero base
	pathLen     uint32               // length of below path
	path        []uint32             // list of nodes in alignment, len of path is = pathLen
	numCigarOps uint32               // number of cigar operations
	byteCigar   []cigar.ByteCigar    // run-length encoding of alignment. Can be copied directly from Giraf struct IF the aligner encodes mismatches into the cigar.
	fancySeq    dnaThreeBit.ThreeBit // three bit encoding of 'fancy-cigar' seq determined at time of compression. Length will be determined by Op's in byteCigar
	alnScore    int64                // alignment quality score
	mapQ        uint8                // mapping quality score
	// Seq field is dropped, can be determined from aln
	numQualOps uint16            // number of qual operations for field below
	qual       []cigar.ByteCigar // phred-scaled base qualities. run-length encoding using ByteCigar
	notes      []binNote         // the notes field will be identical to the SAM notes field. See section 4.2.4 in SAM specs for details
}

type binNote struct {
	tag     [2]byte // tag ID
	tagType byte    // type of data encoded. See section 4.2.4 in SAM specs
	data    string  // data attributed to the the tag. See section 4.2.4 in SAM specs for how to decode
}
