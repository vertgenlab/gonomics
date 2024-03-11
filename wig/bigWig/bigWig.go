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

// specs for version 4 are found in this file: https://github.com/ucscGenomeBrowser/kent/blob/master/src/inc/bbiFile.h

// this is the magic number for bigWig files
const bigWigMagic = 2291137574         // little endian
const bigWigMagicBigEndian = 654086024 // big endian, not supported by this package

// Header represents the bbi header of the bigWig data file and contains metadata and offset values for file navigation.
type Header struct {
	Magic                uint32 // magic bytes at the beginning of file. Must match constant above to be a valid bigWig.
	Version              uint16 // Specifies the wig version. This package was built following specs for bigWig version 4.
	ZoomLevels           uint16 // Number of zoom levels built into the file.
	ChromosomeTreeOffset uint64 // Offset in file to chromosome B+ tree index.
	FullDataOffset       uint64 // Offset to the main dta. Points to the dataCount.
	FullIndexOffset      uint64 // Offset to R tree index of items.
	FieldCount           uint16 // Number of fields in a BED file. For bigWig, this should always be zero.
	DefinedFieldCount    uint16 // Number of fields that are predefined BED files. this should also be zero for bigWig files.
	AutoSqlOffset        uint64 // From specs: Offset to zero-terminated string with .as spec. Used for bigBig, not used in bigWig.
	TotalSummaryOffset   uint64 // Offset to overall file summary data block.
	UncompressBufferSize uint32 // Maximum size of decompression buffer needed (nonzero on compressed files).
	ExtensionOffset      uint64 // Offset to header extension 0 if no such extension.
}

// We currently support only reading bbi headers in LittleEndian.
func readHeader(file *fileio.EasyReader) Header {
	var header = Header{}
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
	err = binary.Read(file, binary.LittleEndian, &header.TotalSummaryOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.UncompressBufferSize)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.ExtensionOffset)
	exception.PanicOnErr(err)
	return header
}
