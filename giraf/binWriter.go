package giraf

import (
	"bytes"
	"encoding/binary"
	"github.com/biogo/hts/bgzf"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaThreeBit"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
)

type BinWriter struct {
	bg  *bgzf.Writer
	buf bytes.Buffer
}

func NewWriter(file io.Writer) *BinWriter {
	return &BinWriter{
		bg: bgzf.NewWriter(file, 1), //TODO: Play with different levels of concurrency
	}
}

func (bw *BinWriter) Write(g *Giraf) error {
	bw.buf.Reset() // clear buffer for new write
	var currBuf [8]byte

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
		bw.buf.Write(currBuf[:2]) // byteCigar.RunLen (uint16)
		bw.buf.Write([]byte{val.Op}) // byteCigar.Op (byte)
	}

	// fancySeq (dnaThreeBit.ThreeBit)
	fancySeq := getFancySeq(g.Seq, g.Cigar)
	binary.LittleEndian.PutUint32(currBuf[:4], uint32(fancySeq.Len))
	bw.buf.Write(currBuf[:4]) // fancySeq.Len (uint32)

	for _, val := range fancySeq.Seq {
		binary.LittleEndian.PutUint64(currBuf[:8], val)
		bw.buf.Write(currBuf[:8]) // fancySeq.Seq (uint64)
	}

	// alnScore (int64)
	binary.LittleEndian.PutUint64(currBuf[:8], uint64(g.AlnScore))
	bw.buf.Write(currBuf[:8])

	// mapQ (uint8)
	bw.buf.Write([]byte{g.MapQ})

	// qual ([]cigar.ByteCigar)
	qual := encodeQual(g.Qual)
	binary.LittleEndian.PutUint16(currBuf[:2], uint16(len(qual)))
	bw.buf.Write(currBuf[:2]) // numQualOps (uint16)

	for _, val := range qual {
		binary.LittleEndian.PutUint16(currBuf[:2], val.RunLen)
		bw.buf.Write(currBuf[:2]) // byteCigar.RunLen (uint16)
		bw.buf.Write([]byte{val.Op}) // byteCigar.Op (byte)
	}

	// notes ([]BinNote)
	for _, val := range g.Notes {
		bw.buf.Write(noteToBytes(val)) // notes ([]BinNotes)
	}

	// write all bytes in buffer
	_, err := bw.bg.Write(bw.buf.Bytes())
	return err
}

func WriteBinGiraf(filename string) {
	inputStream := GoReadToChan(filename)
	outfile := fileio.EasyCreate(filename + ".fe")
	writer := NewWriter(outfile)
	var err error

	for giraf := range inputStream {
		err = writer.Write(giraf)
		common.ExitIfError(err)
	}
}

//func CompressGiraf(g Giraf) BinGiraf {
//	var answer BinGiraf
//	answer.qNameLen = uint8(len(g.QName))
//	answer.qName = g.QName
//	answer.flag = g.Flag
//	answer.tStart = uint32(g.Path.TStart)
//	answer.tEnd = uint32(g.Path.TEnd)
//	answer.pathLen = uint16(len(g.Path.Nodes))
//	answer.path = g.Path.Nodes
//	answer.numCigarOps = uint16(len(g.Cigar))
//	answer.byteCigar = g.Cigar
//	answer.fancySeq = getFancySeq(g.Seq, g.Cigar) //TODO compress ThreeBitDna into []byte + len ?
//	answer.alnScore = int64(g.AlnScore)
//	answer.mapQ = g.MapQ
//	answer.qual = encodeQual(g.Qual)
//	answer.numQualOps = uint16(len(answer.qual))
//	answer.notes = encodeNotes(g.Notes)
//	//TODO hasNamePrefix
//	//TODO blocksize
//	return answer
//}

func getFancySeq(seq []dna.Base, cigar []cigar.ByteCigar) dnaThreeBit.ThreeBit {
	var answer []dna.Base
	var seqIdx int
	for _, val := range cigar {
		if val.Op == 'S' || val.Op == 'X' || val.Op == 'I' {
			answer = append(answer, seq[seqIdx:seqIdx+int(val.RunLen)]...)
		}
		seqIdx += int(val.RunLen)
	}
	return *dnaThreeBit.NewThreeBit(answer, dnaThreeBit.A)
}

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

func noteToBytes(n Note) []byte {
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
		if answer[len(answer)-1] != '\000' { // check if last byte is nil string
			answer = append(answer, '\000') // if not, create nil terminated string
			log.Printf("WARNING: String in tag is not nil terminated, adding nil: %+v", n)
		}
	case 'H': // hex
		if answer[len(answer)-1] != '\000' { // check if last byte is nil string
			answer = append(answer, '\000') // if not, create nil terminated string
			log.Printf("WARNING: String in tag is not nil terminated, adding nil: %+v", n)
		}
	case 'B': // array //TODO: currently treat as string, but could be array per sam format. Will require changes to Giraf struct
		if answer[len(answer)-1] != '\000' { // check if last byte is nil string
			answer = append(answer, '\000') // if not, create nil terminated string
		}
	default:
		log.Fatalf("ERROR: Unrecognized tag type: %s", string(n.Type))
	}
	return answer
}