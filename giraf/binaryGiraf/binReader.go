package binaryGiraf

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"github.com/vertgenlab/gonomics/bgzf"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaThreeBit"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/genomeGraph"
	"github.com/vertgenlab/gonomics/giraf"
	"io"
	"log"
	"strings"
)

// The BinReader struct wraps the bgzf reader from the biogo repository with a bytes buffer to store encoded giraf records
type BinReader struct {
	bg       *bgzf.Reader
	currData *bytes.Buffer
	readBffr []byte
}

// NewBinReader creates a new BinReader
func NewBinReader(filename string) *BinReader {
	return &BinReader{
		bg:       bgzf.NewReader(filename),
		currData: &bytes.Buffer{},
		readBffr: make([]byte, 0),
	}
}

// DecompressGiraf will decode a binary giraf file (.giraf.fe) and output a giraf file (.giraf)
func DecompressGiraf(infilename string, outfilename string, graph *genomeGraph.GenomeGraph) {
	// Initialize infile
	infile := fileio.EasyOpen(infilename)
	defer infile.Close()
	reader := NewBinReader(infilename)
	var err error

	// Initialize outfile
	outfile := fileio.EasyCreate(outfilename)
	defer outfile.Close()

	// Read info until EOF
	var curr giraf.Giraf
	for curr, err = ReadGiraf(reader, graph); err != io.EOF; curr, err = ReadGiraf(reader, graph) {
		common.ExitIfError(err)
		giraf.WriteGirafToFileHandle(outfile, &curr)
	}

	// Close reader
	err = reader.bg.Close()
	common.ExitIfError(err)
}

// The Read method for the BinWriter struct decompresses a single giraf record and writes to file
func ReadGiraf(br *BinReader, g *genomeGraph.GenomeGraph) (giraf.Giraf, error) {
	var answer giraf.Giraf
	var bytesRead int
	var err error
	var i, j, runLen uint16
	var buffer [4]byte

	// blockSize (uint32)
	bytesRead, err = br.bg.Read(buffer[:4])
	if err == io.EOF {
		return answer, io.EOF
	}
	if bytesRead != 4 {
		return answer, io.ErrUnexpectedEOF
	}
	common.ExitIfError(err)
	blockSize := int(binary.LittleEndian.Uint32(buffer[:4]))

	// check cap of readBffr. Either make new slice, or re-slice as needed
	if cap(br.readBffr) < blockSize {
		br.readBffr = make([]byte, blockSize)
	} else {
		br.readBffr = br.readBffr[:blockSize]
	}

	// read block to readBffr and fill currData buffer with read bytes
	bytesRead, err = io.ReadFull(br.bg, br.readBffr)
	if err == io.ErrUnexpectedEOF || bytesRead != blockSize {
		log.Fatalln("ERROR: truncated record")
	}
	common.ExitIfError(err)
	br.currData.Write(br.readBffr) // copy data from intermediate buffer to bytes buffer

	// qNameLen (uint8)
	qNameLen := int(br.currData.Next(1)[0])

	// qName (string)
	answer.QName = string(br.currData.Next(qNameLen))

	// flag (uint8)
	answer.Flag = br.currData.Next(1)[0]

	// tStart (uint32)
	answer.Path.TStart = int(binary.LittleEndian.Uint32(br.currData.Next(4)))

	// tEnd (uint32)
	answer.Path.TEnd = int(binary.LittleEndian.Uint32(br.currData.Next(4)))

	// path ([]uint32
	pathLen := binary.LittleEndian.Uint32(br.currData.Next(4)) // pathLen (uint32)
	var node, k uint32
	for k = 0; k < pathLen; k++ { // path ([]uint32)
		node = binary.LittleEndian.Uint32(br.currData.Next(4))
		answer.Path.Nodes = append(answer.Path.Nodes, node)
	}

	// cigar ([]cigar.ByteCigar)
	numCigarOps := binary.LittleEndian.Uint32(br.currData.Next(4)) // numCigarOps (uint16)
	for k = 0; k < numCigarOps; k++ {
		answer.Cigar = append(answer.Cigar, cigar.ByteCigar{
			RunLen: binary.LittleEndian.Uint16(br.currData.Next(2)), // byteCigar.RunLen (uint16)
			Op:     br.currData.Next(1)[0],                          // byteCigar.Op (byte)
		})
	}

	// fancySeq (dnaThreeBit.ThreeBit)
	var fancyBases dnaThreeBit.ThreeBit
	fancyBases.Len = int(binary.LittleEndian.Uint32(br.currData.Next(4))) // fancySeq.Len (uint32)

	// figure out how many uint64 are needed to store fancySeqLen worth of bases
	numInts := fancyBases.Len / 21
	if fancyBases.Len%21 != 0 {
		numInts++
	}

	for j := 0; j < numInts; j++ {
		fancyBases.Seq = append(fancyBases.Seq, binary.LittleEndian.Uint64(br.currData.Next(8))) // fancySeq.Seq (uint64)
	}
	addFullSeq(&answer, &fancyBases, g)

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

	// QStart, QEnd (int)
	answer.QStart, answer.QEnd = determineQStartQEnd(&answer)

	// notes ([]BinNote)
	appendNotes(&answer, br) // notes ([]bytes)

	if giraf.IsForwardRead(&answer) {
		answer.PosStrand = true
	}

	return answer, err
}

// addFullSeq parses the cigar and the fancySeq fields to retrieve the full length read sequence
func addFullSeq(answer *giraf.Giraf, fancySeq *dnaThreeBit.ThreeBit, graph *genomeGraph.GenomeGraph) {
	var fancyBases []dna.Base
	if fancySeq.Len != 0 {
		fancyBases = dnaThreeBit.ToDnaBases(fancySeq)
	}
	refIdx := answer.Path.TStart
	var currNodeId int
	var currNode genomeGraph.Node = graph.Nodes[answer.Path.Nodes[0]]
	var i uint16

	for _, cigar := range answer.Cigar {
		switch cigar.Op {
		case '=':
			for i = 0; i < cigar.RunLen; i++ {
				if refIdx > len(currNode.Seq)-1 {
					refIdx = 0
					currNodeId++
					currNode = graph.Nodes[answer.Path.Nodes[currNodeId]]
				}

				answer.Seq = append(answer.Seq, currNode.Seq[refIdx])
				refIdx++
			}
		case 'X':
			answer.Seq = append(answer.Seq, fancyBases[:cigar.RunLen]...) // retrieve RunLen number of bases from the fancyBases slice
			fancyBases = fancyBases[:cigar.RunLen]                        // removed the retrieved bases from the fancyBases slice
			refIdx += int(cigar.RunLen)
		case 'S':
			answer.Seq = append(answer.Seq, fancyBases[:cigar.RunLen]...) // retrieve RunLen number of bases from the fancyBases slice
			fancyBases = fancyBases[:cigar.RunLen]                        // removed the retrieved bases from the fancyBases slice
		case 'I':
			answer.Seq = append(answer.Seq, fancyBases[:cigar.RunLen]...) // retrieve RunLen number of bases from the fancyBases slice
			fancyBases = fancyBases[:cigar.RunLen]                        // removed the retrieved bases from the fancyBases slice
		case 'D':
			refIdx += int(cigar.RunLen)
		default:
			log.Fatalf("ERROR: Unrecognized cigar operation: %v", cigar.Op)
		}
	}

}

// appendNotes parses the encoded notes at the end of the alignment record and appends them to the output giraf record
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

// determineQStartQEnd parses the cigar to determine where the alignment starts and ends
func determineQStartQEnd(answer *giraf.Giraf) (int, int) {
	var start, end int
	if answer.Cigar == nil {
		return 0, 0
	}

	if answer.Cigar[0].Op == 'S' {
		start = int(answer.Cigar[0].RunLen)
	} else {
		start = 0
	}

	if answer.Cigar[len(answer.Cigar)-1].Op == 'S' {
		end = (len(answer.Seq) - 1) - int(answer.Cigar[len(answer.Cigar)-1].RunLen)
	} else {
		end = len(answer.Seq) - 1
	}

	return start, end
}
