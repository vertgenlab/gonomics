package binaryGiraf

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"github.com/biogo/hts/bgzf"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"io"
	"log"
	"strings"
)

type BinReader struct {
	bg  *bgzf.Reader
	currData *bytes.Buffer
}

func NewBinReader(file io.Reader) *BinReader {
	reader, err := bgzf.NewReader(file, 1) //TODO: Play with different levels of concurrency
	common.ExitIfError(err)
	return &BinReader{
		bg: reader,
		currData: &bytes.Buffer{},
	}
}

func DecompressGiraf(filename string, graph *simpleGraph.SimpleGraph) {
	// Initialize infile
	infile := fileio.EasyOpen(filename)
	defer infile.Close()
	reader := NewBinReader(infile.BuffReader)
	var err error

	// Initialize outfile
	outfile := fileio.EasyCreate(strings.TrimSuffix(filename, ".fe"))
	defer outfile.Close()

	// Read info until EOF
	var curr giraf.Giraf
	for curr, err = reader.Read(graph); err != io.EOF; curr, err = reader.Read(graph) {
		common.ExitIfError(err)
		giraf.WriteGirafToFileHandle(outfile, &curr)
	}

	// Close reader
	err = reader.bg.Close()
	common.ExitIfError(err)
}

func (br *BinReader) Read(g *simpleGraph.SimpleGraph) (giraf.Giraf, error) {
	var answer giraf.Giraf
	var bytesRead int
	var err error
	buf := tempBuf{
		reader: br,
	}

	// reset data buffer
	br.currData.Reset()

	// blockSize (uint32)
	bytesRead, err = buf.reader.bg.Read(buf.buffer[:4])
	if err == io.EOF || bytesRead != 4 {
		return answer, io.EOF
	}
	common.ExitIfError(err)
	blockSize := int(binary.LittleEndian.Uint32(buf.buffer[:4]))

	// fill currData buffer with blockSize of bytes
	data := make([]byte, blockSize)
	bytesRead, err = io.ReadFull(br.bg, data)
	if err == io.EOF || bytesRead != blockSize {
		log.Fatalln("ERROR: truncated record")
	}
	common.ExitIfError(err)
	br.currData.Write(data)

	fmt.Println(br.currData.Bytes())

	// qNameLen (uint8)

	// qName (string)

	// flag (uint8)

	// tStart (uint32)

	// tEnd (uint32)

	// path ([]uint32
		// pathLen (uint16)
		// path ([]uint32)

	// cigar ([]cigar.ByteCigar)
		// numCigarOps (uint16)
		// byteCigar.RunLen (uint16)
		// byteCigar.Op (byte)

	// fancySeq (dnaThreeBit.ThreeBit)
		// fancySeq.Len (uint32)
		// fancySeq.Seq (uint64)

	// alnScore (int64)

	// mapQ (uint8)

	// qual ([]cigar.ByteCigar)
		// numQualOps (uint16)
		// byteCigar.RunLen (uint16)
		// byteCigar.Op (byte)

	// notes ([]BinNote)
		// notes ([]bytes)

	return answer, err
}

type tempBuf struct {
	buffer  [8]byte
	reader	*BinReader
}

func (b tempBuf) readUint8() uint8 {
	bytesRead, err := b.reader.bg.Read(b.buffer[:1])
	if err == io.EOF || bytesRead != 1 {
		log.Fatalln("ERROR: truncated record")
	}
	common.ExitIfError(err)
	return b.buffer[0]
}

func (b tempBuf) readUint16() uint16 {
	bytesRead, err := b.reader.bg.Read(b.buffer[:2])
	if err == io.EOF || bytesRead != 2 {
		log.Fatalln("ERROR: truncated record")
	}
	common.ExitIfError(err)
	return binary.LittleEndian.Uint16(b.buffer[:2])
}

func (b tempBuf) readUint32() uint32 {
	bytesRead, err := b.reader.bg.Read(b.buffer[:4])
	if err == io.EOF || bytesRead != 4 {
		log.Fatalln("ERROR: truncated record")
	}
	common.ExitIfError(err)
	return binary.LittleEndian.Uint32(b.buffer[:4])
}

func (b tempBuf) readUint64() uint64 {
	bytesRead, err := b.reader.bg.Read(b.buffer[:8])
	if err == io.EOF || bytesRead != 8 {
		log.Fatalln("ERROR: truncated record")
	}
	common.ExitIfError(err)
	return binary.LittleEndian.Uint64(b.buffer[:8])
}
