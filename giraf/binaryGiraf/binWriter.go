package binaryGiraf

import (
	"bytes"
	"encoding/binary"
	"github.com/biogo/hts/bgzf"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaThreeBit"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"io"
	"log"
)

// The BinWriter struct wraps the bgzf writer from the biogo repository with a bytes buffer to store encoded giraf records
type BinWriter struct {
	bg  *bgzf.Writer
	buf bytes.Buffer
}

// NewBinWriter creates a new BinWriter
func NewBinWriter(file io.Writer) *BinWriter {
	return &BinWriter{
		bg: bgzf.NewWriter(file, 1), //TODO: Play with different levels of concurrency
	}
}

// CompressGiraf will encode a giraf file (.giraf) and output a binary giraf file (.giraf.fe)
func CompressGiraf(filename string) {
	inputStream := giraf.GoReadToChan(filename)
	outfile := fileio.EasyCreate(filename + ".fe")
	defer outfile.Close()
	writer := NewBinWriter(outfile)
	var err error

	// Write info from all girafs in inputStream
	for giraf := range inputStream {
		err = writer.Write(giraf)
		common.ExitIfError(err)
	}

	// Close writer
	err = writer.bg.Close()
	common.ExitIfError(err)
}

// Byte size of BinGiraf fixed size fields excluding blockSize.
// Exluded Fields: qName, path, byteCigar, fancySeq.Seq, qual, notes
var binGirafFixedSize int = 29

// The Write method for the BinWriter struct compresses a single giraf record and writes to file
func (bw *BinWriter) Write(g *giraf.Giraf) error {
	bw.buf.Reset() // clear buffer for new write
	var currBuf [8]byte

	// encode data for later use
	fancySeq := getFancySeq(g.Seq, g.Cigar)
	qual := encodeQual(g.Qual)
	notes := notesToBytes(g.Notes)

	// blockSize (uint32)
	size := binGirafFixedSize +
		len(g.QName) +
		(4 * len(g.Path.Nodes)) +
		(3 * len(g.Cigar)) +
		(8 * len(fancySeq.Seq)) +
		(3 * len(qual)) +
		len(notes)
	binary.LittleEndian.PutUint32(currBuf[:4], uint32(size))
	bw.buf.Write(currBuf[:4])

	// qNameLen (uint8)
	bw.buf.Write([]byte{uint8(len(g.QName))})

	// qName (string)
	bw.buf.Write([]byte(g.QName))

	// flag (uint8)
	bw.buf.Write([]byte{g.Flag})

	// tStart (uint32)
	binary.LittleEndian.PutUint32(currBuf[:4], uint32(g.Path.TStart))
	bw.buf.Write(currBuf[:4])

	// tEnd (uint32)
	binary.LittleEndian.PutUint32(currBuf[:4], uint32(g.Path.TEnd))
	bw.buf.Write(currBuf[:4])

	// path ([]uint32
	binary.LittleEndian.PutUint16(currBuf[:2], uint16(len(g.Path.Nodes)))
	bw.buf.Write(currBuf[:2]) // pathLen (uint16)

	for _, val := range g.Path.Nodes {
		binary.LittleEndian.PutUint32(currBuf[:4], val)
		bw.buf.Write(currBuf[:4]) // path ([]uint32)
	}

	// cigar ([]cigar.ByteCigar)
	binary.LittleEndian.PutUint16(currBuf[:2], uint16(len(g.Cigar)))
	bw.buf.Write(currBuf[:2]) // numCigarOps (uint16)

	for _, val := range g.Cigar {
		binary.LittleEndian.PutUint16(currBuf[:2], val.RunLen)
		bw.buf.Write(currBuf[:2])    // byteCigar.RunLen (uint16)
		bw.buf.Write([]byte{val.Op}) // byteCigar.Op (byte)
	}

	// fancySeq (dnaThreeBit.ThreeBit)
	binary.LittleEndian.PutUint32(currBuf[:4], uint32(fancySeq.Len))
	bw.buf.Write(currBuf[:4]) // fancySeq.Len (uint32)

	for _, val := range fancySeq.Seq {
		binary.LittleEndian.PutUint64(currBuf[:8], val) // TODO: does not need to be uint64?
		bw.buf.Write(currBuf[:8])                       // fancySeq.Seq (uint64)
	}

	// alnScore (int64)
	binary.LittleEndian.PutUint64(currBuf[:8], uint64(g.AlnScore))
	bw.buf.Write(currBuf[:8])

	// mapQ (uint8)
	bw.buf.Write([]byte{g.MapQ})

	// qual ([]cigar.ByteCigar)
	binary.LittleEndian.PutUint16(currBuf[:2], uint16(len(qual)))
	bw.buf.Write(currBuf[:2]) // numQualOps (uint16)

	for _, val := range qual {
		binary.LittleEndian.PutUint16(currBuf[:2], val.RunLen)
		bw.buf.Write(currBuf[:2])    // byteCigar.RunLen (uint16)
		bw.buf.Write([]byte{val.Op}) // byteCigar.Op (byte)
	}

	// notes ([]BinNote)
	bw.buf.Write(notes) // notes ([]bytes)

	// write all bytes in buffer
	_, err := bw.bg.Write(bw.buf.Bytes())
	return err
}

// getFancySeq will parse the []cigar.ByteCigar and record any bases that cannot be recovered by reference matching
func getFancySeq(seq []dna.Base, cigar []cigar.ByteCigar) dnaThreeBit.ThreeBit {
	var answer []dna.Base
	var seqIdx int
	if cigar == nil {
		return *dnaThreeBit.NewThreeBit(seq, dnaThreeBit.A)
	}
	for _, val := range cigar {
		if val.Op == 'S' || val.Op == 'X' || val.Op == 'I' {
			answer = append(answer, seq[seqIdx:seqIdx+int(val.RunLen)]...)
		}
		seqIdx += int(val.RunLen)
	}
	return *dnaThreeBit.NewThreeBit(answer, dnaThreeBit.A)
}

// encodeQual creates a run-length encoded representation of the Qual scores in the ByteCigar format
func encodeQual(q []uint8) []cigar.ByteCigar {
	answer := make([]cigar.ByteCigar, 0, len(q))
	var curr cigar.ByteCigar
	curr.Op = q[0]
	for i := 0; i < len(q); i++ {
		if q[i] != curr.Op && curr.RunLen != 0 {
			answer = append(answer, curr)
			curr.RunLen = 0
			curr.Op = q[i]
		}
		curr.RunLen++
	}

	if curr.RunLen != 0 {
		answer = append(answer, curr)
	}

	return answer
}

// notesToBytes parses a slice of notes and returns a slice of encoded bytes
func notesToBytes(n []giraf.Note) []byte {
	var answer []byte
	for _, val := range n {
		answer = append(answer, noteToBytes(val)...)
	}
	return answer
}

// noteToBytes parses a single note and returns a slice of encoded bytes
func noteToBytes(n giraf.Note) []byte {
	var answer []byte
	var currBuf [4]byte
	if len(n.Tag) != 2 {
		log.Fatalf("ERROR: Tag is not two bytes: %s", n.Tag)
	}
	answer = append(answer, n.Tag[0], n.Tag[1], n.Type)
	switch n.Type {
	case 'A': // rune
		answer = append(answer, n.Value...)
		if len(answer) != 4 {
			log.Fatalf("ERROR: Improperly formatted note: %+v", n)
		}
	case 'c': // int8
		val := common.StringToInt8(n.Value)
		answer = append(answer, byte(val))
		if len(answer) != 4 {
			log.Fatalf("ERROR: Improperly formatted note: %+v", n)
		}
	case 'C': // uint8
		val := common.StringToUint8(n.Value)
		answer = append(answer, byte(val))
		if len(answer) != 4 {
			log.Fatalf("ERROR: Improperly formatted note: %+v", n)
		}
	case 's': // int16
		val := common.StringToInt16(n.Value)
		binary.LittleEndian.PutUint16(currBuf[:2], uint16(val))
		answer = append(answer, currBuf[:2]...)
		if len(answer) != 5 {
			log.Fatalf("ERROR: Improperly formatted note: %+v", n)
		}
	case 'S': // uint16
		val := common.StringToUint16(n.Value)
		binary.LittleEndian.PutUint16(currBuf[:2], val)
		answer = append(answer, currBuf[:2]...)
		if len(answer) != 5 {
			log.Fatalf("ERROR: Improperly formatted note: %+v", n)
		}
	case 'i': // int32
		val := common.StringToInt32(n.Value)
		binary.LittleEndian.PutUint32(currBuf[:4], uint32(val))
		answer = append(answer, currBuf[:4]...)
		if len(answer) != 7 {
			log.Fatalf("ERROR: Improperly formatted note: %+v", n)
		}
	case 'I': // uint32
		val := common.StringToUint32(n.Value)
		binary.LittleEndian.PutUint32(currBuf[:4], val)
		answer = append(answer, currBuf[:4]...)
		if len(answer) != 7 {
			log.Fatalf("ERROR: Improperly formatted note: %+v", n)
		}
	case 'f': // float32
		val := uint32(common.StringToFloat64(n.Value))
		binary.LittleEndian.PutUint32(currBuf[:4], val)
		answer = append(answer, currBuf[:4]...)
		if len(answer) != 7 {
			log.Fatalf("ERROR: Improperly formatted note: %+v", n)
		}
	case 'Z': // string
		answer = append(answer, n.Value...)
		if answer[len(answer)-1] != '\000' { // check if last byte is nil string
			answer = append(answer, '\000') // if not, create nil terminated string
		}
	case 'H': // hex
		answer = append(answer, n.Value...)
		if answer[len(answer)-1] != '\000' { // check if last byte is nil string
			answer = append(answer, '\000') // if not, create nil terminated string
		}
	case 'B': // array //TODO: currently treat as string, but could be array per sam format. Will require changes to Giraf struct
		answer = append(answer, n.Value...)
		if answer[len(answer)-1] != '\000' { // check if last byte is nil string
			answer = append(answer, '\000') // if not, create nil terminated string
		}
	default:
		log.Fatalf("ERROR: Unrecognized tag type: %s", string(n.Type))
	}
	return answer
}
