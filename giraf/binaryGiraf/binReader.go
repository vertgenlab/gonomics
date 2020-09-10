package binaryGiraf

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"github.com/biogo/hts/bgzf"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaThreeBit"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"io"
	"log"
	"strings"
)

type BinReader struct {
	bg       *bgzf.Reader
	currData *bytes.Buffer
}

func NewBinReader(file io.Reader) *BinReader {
	reader, err := bgzf.NewReader(file, 1) //TODO: Play with different levels of concurrency
	common.ExitIfError(err)
	return &BinReader{
		bg:       reader,
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
	var i, j, runLen uint16
	var buffer [4]byte

	// reset data buffer
	br.currData.Reset()

	// blockSize (uint32)
	bytesRead, err = br.bg.Read(buffer[:4])
	if err == io.EOF || bytesRead != 4 {
		return answer, io.EOF
	}
	common.ExitIfError(err)
	blockSize := int(binary.LittleEndian.Uint32(buffer[:4]))

	// fill currData buffer with blockSize of bytes
	data := make([]byte, blockSize) // read giraf to intermediate buffer //TODO: this is not great, find a way to not make an intermediate buffer
	bytesRead, err = io.ReadFull(br.bg, data)
	if err == io.EOF || bytesRead != blockSize {
		log.Fatalln("ERROR: truncated record")
	}
	common.ExitIfError(err)
	br.currData.Write(data) // copy data from intermediate buffer to bytes buffer

	// qNameLen (uint8)
	qNameLen := int(br.currData.Next(1)[0])

	// qName (string)
	answer.QName = string(br.currData.Next(qNameLen))

	// flag (uint8)
	answer.Flag = br.currData.Next(1)[0]

	// tStart (uint32)
	answer.QStart = int(binary.LittleEndian.Uint32(br.currData.Next(4)))
	answer.Path.TStart = answer.QStart //TODO: Is this true?

	// tEnd (uint32)
	answer.QEnd = int(binary.LittleEndian.Uint32(br.currData.Next(4)))
	answer.Path.TEnd = answer.QEnd //TODO: Is this true?

	// path ([]uint32
	pathLen := binary.LittleEndian.Uint16(br.currData.Next(2)) // pathLen (uint16)
	var node uint32
	for i = 0; i < pathLen; i++ { // path ([]uint32)
		node = binary.LittleEndian.Uint32(br.currData.Next(4))
		answer.Path.Nodes = append(answer.Path.Nodes, node)
	}

	// cigar ([]cigar.ByteCigar)
	numCigarOps := binary.LittleEndian.Uint16(br.currData.Next(2)) // numCigarOps (uint16)
	for i = 0; i < numCigarOps; i++ {
		answer.Cigar = append(answer.Cigar, cigar.ByteCigar{
			RunLen: binary.LittleEndian.Uint16(br.currData.Next(2)), // byteCigar.RunLen (uint16)
			Op:     br.currData.Next(1)[0],                          // byteCigar.Op (byte)
		})
	}

	// fancySeq (dnaThreeBit.ThreeBit)
	var fancyBases dnaThreeBit.ThreeBit
	fancyBases.Len = int(binary.LittleEndian.Uint32(br.currData.Next(4))) // fancySeq.Len (uint32)

	// figure out how many uint64 are needed to store fancySeqLen worth of bases
	//TODO: Craig, do you think this is the best way to determine this?
	numInts := fancyBases.Len / 21
	if fancyBases.Len%21 != 0 {
		numInts++
	}

	for j := 0; j < numInts; j++ {
		fancyBases.Seq = append(fancyBases.Seq, binary.LittleEndian.Uint64(br.currData.Next(8))) // fancySeq.Seq (uint64)
	}

	answer.Seq = fancyToFullSeq(&fancyBases, g)

	// alnScore (int64)
	answer.AlnScore = int(binary.LittleEndian.Uint64(br.currData.Next(8)))

	// mapQ (uint8)
	answer.MapQ = br.currData.Next(1)[0]

	// qual ([]cigar.ByteCigar)
	var op uint8
	numQualOps := binary.LittleEndian.Uint16(br.currData.Next(2)) // numQualOps (uint16)
	for i = 0; i < numQualOps; i++ {
		runLen = binary.LittleEndian.Uint16(br.currData.Next(2)) // byteCigar.RunLen (uint16)
		op = br.currData.Next(1)[0]                              // byteCigar.Op (byte)
		for j = 0; j < runLen; j++ {
			answer.Qual = append(answer.Qual, op)
		}
	}

	// notes ([]BinNote)
	appendNotes(&answer, br) // notes ([]bytes)

	return answer, err
}

func fancyToFullSeq(fancySeq *dnaThreeBit.ThreeBit, graph *simpleGraph.SimpleGraph) []dna.Base {
	answer := make([]dna.Base, 0)
	//TODO
	return answer
}

func appendNotes(answer *giraf.Giraf, br *BinReader) {
	var currNote giraf.Note
	var currString strings.Builder
	for br.currData.Len() != 0 {
		currNote.Tag = br.currData.Next(2)
		currNote.Type = br.currData.Next(1)[0]
		switch currNote.Type {
		case 'A': // rune
			currNote.Value = string(rune(br.currData.Next(1)[0]))

		case 'c': // int8
			currNote.Value = fmt.Sprintf("%d", int8(br.currData.Next(1)[0]))

		case 'C': // uint8
			currNote.Value = fmt.Sprintf("%d", br.currData.Next(1)[0])

		case 's': // int16
			currNote.Value = fmt.Sprintf("%d", int16(binary.LittleEndian.Uint16(br.currData.Next(2))))

		case 'S': // uint16
			currNote.Value = fmt.Sprintf("%d", binary.LittleEndian.Uint16(br.currData.Next(2)))

		case 'i': // int32
			currNote.Value = fmt.Sprintf("%d", int32(binary.LittleEndian.Uint32(br.currData.Next(4))))

		case 'I': // uint32
			currNote.Value = fmt.Sprintf("%d", binary.LittleEndian.Uint32(br.currData.Next(4)))

		case 'f': // float32
			currNote.Value = fmt.Sprintf("%f", float32(binary.LittleEndian.Uint32(br.currData.Next(4))))

		case 'Z': // string
			currString.Reset()
			currString.WriteByte(br.currData.Next(1)[0])
			for currString.String()[currString.Len()-1] != '\000' {
				currString.WriteByte(br.currData.Next(1)[0])
			}
			currNote.Value = currString.String()[:currString.Len()-1]

		case 'H': // hex
			currString.Reset()
			currString.WriteByte(br.currData.Next(1)[0])
			for currString.String()[currString.Len()-1] != '\000' {
				currString.WriteByte(br.currData.Next(1)[0])
			}
			currNote.Value = currString.String()[:currString.Len()-1]

		case 'B': // array //TODO: currently treat as string, but could be array per sam format. Will require changes to Giraf struct
			currString.Reset()
			currString.WriteByte(br.currData.Next(1)[0])
			for currString.String()[currString.Len()-1] != '\000' {
				currString.WriteByte(br.currData.Next(1)[0])
			}
			currNote.Value = currString.String()[:currString.Len()-1]

		default:
			log.Fatalf("ERROR: Unrecognized tag type: %s", string(currNote.Type))
		}
		answer.Notes = append(answer.Notes, currNote)
	}
}
