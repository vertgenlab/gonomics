package bigWig

import (
	"encoding/binary"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

// the official documentation for the bigWig file format
// can be found in the supplementary materials of the following
// publication:
// Kent, W.J., et al. BigWig and BigBed: enabling browsing of large distributed
//datasets. (2010) Bioinformatics, 26(17) 2204-2207.

// this is the magic number for bigWig files
const bigWigMagic = 2291137574

type Header struct {
	Magic                uint32 // magic bytes at the beginning of file. Must match constant above to be a valid bigWig.
	Version              uint16 // Specifies the wig version. This package was built following specs for bigWig version 3.
	ZoomLevels           uint16 // Number of zoom levels built into the file.
	ChromosomeTreeOffset uint64 // Offset in file to chromosome B+ tree index.
	FullDataOffset       uint64 // Offset to the main dta. Points to the dataCount.
	FullIndexOffset      uint64 // Offset to R tree index of items.
	FieldCount           uint16 // Number of fields in a BED file. For bigWig, this should always be zero.
	DefinedFieldCount    uint16 // Number of fields that are predefined BED files. this should also be zero for bigWig files.
	AutoSqlOffset        uint64 // From specs: Offset to zero-terminated string with .as spec. Used for bigBig, not used in bigWig.
	TotalSummaryOffset   uint64 // Offset to overall file summary data block.
	UncompressBufferSize uint32 // Maximum size of decompression buffer needed (nonzero on compressed files).
	Reserved             uint64 //Unused values, reserved for future expansion. Should always be zero in a valid bigWig.
}

func readHeader(file *fileio.EasyReader) Header {
	var header = Header{}
	err := binary.Read(file, binary.LittleEndian, &header.Magic)
	exception.PanicOnErr(err)
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
	err = binary.Read(file, binary.LittleEndian, &header.TotalSummaryOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.UncompressBufferSize)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.Reserved)
	exception.PanicOnErr(err)
	if header.Reserved != 0 {
		log.Fatalf("Error: bigWig header reserved field must be zero. Found: %v.\n", header.Reserved)
	}
	return header
}
