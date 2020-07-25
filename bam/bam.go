package bam

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"io"
	"log"
	"os"
	"strings"
)

//Note: This is not a bam struct but palce holders for pointers to decode binary file
//4.2 The BAM format defined in hts specs

type BamData struct {
	RName     int32
	Pos       int32
	Bai       uint16
	MapQ      uint8
	RNLength  uint8
	Flag      uint16
	NCigarOp  uint16
	LSeq      int32
	NextRefID int32
	NextPos   int32
	TLength   int32
	QName     string
	Cigar     []uint32
	Seq       []byte
	Qual      []byte
	Aux       []*BamAuxiliary
}

type BamHeader struct {
	Txt       string
	ChromSize map[string]int
	decoder   *cInfo
	Chroms    []*chromInfo.ChromInfo
}

type BamReader struct {
	File   *os.File
	gunzip *Bgzip
	Data      []byte
	BytesRead int
	Error     error
}


type BamAuxiliary struct {
	Tag   [2]byte
	Type  byte
	Value interface{}
}

func BamBlockToSam(header *BamHeader, bam *BamData) *sam.SamAln {
	return &sam.SamAln{
		QName: bam.QName,
		Flag:  int64(bam.Flag),
		RName: header.Chroms[bam.RName].Name,
		Pos:   int64(bam.Pos) + 1,
		MapQ:  int64(bam.MapQ),
		Cigar: ToHeavyCigar(bam.Cigar),
		RNext: setRNext(header, bam),
		PNext: int64(bam.NextPos) + 1,
		TLen:  int64(bam.TLength),
		Seq:   dna.StringToBases(bamSequence(bam.Seq)),
		Qual:  qualToString(bam.Qual),
		Extra: NotesToString(bam.Aux),
	}
}

func setRNext(header *BamHeader, bam *BamData) string {
	if bam.NextRefID == bam.RName {
		return "="
	} else if bam.NextRefID > 0 {
		return header.Chroms[bam.NextRefID].Name //fmt.Sprintf("Chrom: %d", bam.NextRefID)
	} else {
		return "*"
	}
}

func bamSequence(seq []byte) string {
	var buffer bytes.Buffer
	writer := bufio.NewWriter(&buffer)
	t := []byte{'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}
	for i := 0; i < len(seq); i++ {
		b1 := seq[i] >> 4
		b2 := seq[i] & 0xf
		fmt.Fprintf(writer, "%c", t[b1])
		if b2 != 0 {
			fmt.Fprintf(writer, "%c", t[b2])
		}
	}
	writer.Flush()
	return buffer.String()
}

func qualToString(qual []byte) string {
	var buffer bytes.Buffer
	writer := bufio.NewWriter(&buffer)
	fmt.Fprintf(writer, "%s", string(formatQual(qual)))
	writer.Flush()
	return buffer.String()
}

func formatQual(q []byte) []byte {
	for _, v := range q {
		if v != 0xff {
			a := make([]byte, len(q))
			for i, p := range q {
				a[i] = p + 33
			}
			return a
		}
	}
	return []byte{'*'}
}

func NewBamReader(filename string) *BamReader {
	var bamR *BamReader = &BamReader{}
	var err error
	bamR.File = fileio.MustOpen(filename)
	bamR.gunzip = NewBgzipReader(bamR.File)
	bamR.Error = err
	return bamR
}

func ReadHeader(reader *BamReader) *BamHeader {
	bamHeader := MakeHeader()
	magic := make([]byte, 4)
	magic = BgzipBuffer(reader.gunzip, magic)
	if string(magic) != "BAM\001" {
		log.Fatalf("Not a BAM file: %s\n", string(reader.Data))
	}
	reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &bamHeader.decoder.Len)
	common.ExitIfError(reader.Error)
	reader.Data = make([]byte, bamHeader.decoder.Len)
	var i, j, k int = 0, 0, 0
	for i = 0; i < int(bamHeader.decoder.Len); {
		j, reader.Error = reader.gunzip.Read(reader.Data[i:])
		common.ExitIfError(reader.Error)
		i += j
	}
	bamHeader.Txt = string(reader.Data)
	reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &bamHeader.decoder.NRef)
	common.ExitIfError(reader.Error)
	reader.Data = make([]byte, bamHeader.decoder.NRef)
	var lengthName, lengthSeq int32
	for i = 0; i < int(bamHeader.decoder.NRef); i++ {
		lengthName, lengthSeq = 0, 0
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &lengthName)
		common.ExitIfError(reader.Error)
		reader.Data = make([]byte, lengthName)
		for j = 0; j < int(lengthName); {
			k, reader.Error = reader.gunzip.Read(reader.Data[j:])
			common.ExitIfError(reader.Error)
			j += k
		}
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &lengthSeq)
		common.ExitIfError(reader.Error)
		bamHeader.Chroms = append(bamHeader.Chroms, &chromInfo.ChromInfo{Name: strings.Trim(string(reader.Data), "\n\000"), Size: int64(lengthSeq), Order: int64(len(bamHeader.Chroms))})
	}
	return bamHeader
}
/*
func decodeBinaryInt32(reader *BamReader) int {
	var key int32 = 0
	reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &key)
	reader.Data = make([]byte, key)
	var i, j int = 0, 0
	for i = 0; i < int(key); {
		j, reader.Error = reader.gunzip.Read(reader.Data[i:])
		common.ExitIfError(reader.Error)
		i += j
	}
	return i
}*/

func Read(filename string) []*sam.SamAln {
	bamFile := NewBamReader(filename)
	defer bamFile.File.Close()
	h := ReadHeader(bamFile)
	binaryData := make(chan *BamData)
	var ans []*sam.SamAln
	go bamToChan(bamFile, binaryData)
	for each := range binaryData {
		ans = append(ans, BamBlockToSam(h, each))
	}
	return ans
}

func bamToChan(reader *BamReader, binaryData chan<- *BamData) {
	var blockSize int32
	var bitFlag uint32
	var stats uint32
	var i, j int
	var b byte
	var block *BamData

	for {
		block = &BamData{}
		buf := bytes.NewBuffer([]byte{})
		// read block size
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &blockSize)
		if reader.Error == io.EOF {
			close(binaryData)
			break
		}
		common.ExitIfError(reader.Error)
		//read block data
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &block.RName)
		common.ExitIfError(reader.Error)

		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &block.Pos)
		common.ExitIfError(reader.Error)

		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &stats)
		common.ExitIfError(reader.Error)

		block.Bai = uint16((stats >> 16) & 0xffff)
		block.MapQ = uint8((stats >> 8) & 0xff)
		block.RNLength = uint8((stats >> 0) & 0xff)

		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &bitFlag)
		common.ExitIfError(reader.Error)

		// get Flag and NCigarOp from bitFlag
		block.Flag = uint16(bitFlag >> 16)
		block.NCigarOp = uint16(bitFlag & 0xffff)

		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &block.LSeq)
		common.ExitIfError(reader.Error)

		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextRefID)
		common.ExitIfError(reader.Error)

		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextPos)
		common.ExitIfError(reader.Error)

		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &block.TLength)
		common.ExitIfError(reader.Error)

		// parse the read name

		for {
			reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &b)
			common.ExitIfError(reader.Error)

			if b == 0 {
				block.QName = buf.String()
				break
			}
			buf.WriteByte(b)
		}
		// parse cigar block
		block.Cigar = make([]uint32, block.NCigarOp)
		for i = 0; i < int(block.NCigarOp) && reader.Error == nil; i++ {
			reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &block.Cigar[i])
			common.ExitIfError(reader.Error)

		}
		// parse seq
		block.Seq = make([]byte, (block.LSeq+1)/2)
		for i = 0; i < int((block.LSeq+1)/2); i++ {
			reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &block.Seq[i])
			common.ExitIfError(reader.Error)

		}
		block.Qual = make([]byte, block.LSeq)
		for i = 0; i < int(block.LSeq); i++ {
			reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &block.Qual[i])
			common.ExitIfError(reader.Error)
		}
		// read auxiliary data
		j = 8*4 + int(block.RNLength) + 4*int(block.NCigarOp) + int((block.LSeq+1)/2) + int(block.LSeq)
		for i = 0; j+i < int(blockSize); {
			block.Aux = append(block.Aux, ReadBamAuxiliary(reader))
			i += reader.BytesRead
		}
		binaryData <- block
	}
}

func NotesToString(aux []*BamAuxiliary) string {
	var ans []string
	for i := 0; i < len(aux); i++ {
		ans = append(ans, fmt.Sprintf("%c%c:%c:%v", aux[i].Tag[0], aux[i].Tag[1], aux[i].Type, aux[i].Value))
	}
	return strings.Join(ans, "\t")
}

func ReadBamAuxiliary(reader *BamReader) *BamAuxiliary {
	aux := &BamAuxiliary{}
	var i int
	// number of read bytes
	reader.BytesRead = 0
	// read data
	reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &aux.Tag[0])
	common.ExitIfError(reader.Error)

	reader.BytesRead += 1
	reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &aux.Tag[1])
	common.ExitIfError(reader.Error)

	reader.BytesRead += 1
	reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &aux.Type)
	common.ExitIfError(reader.Error)

	// three bytes read so far
	reader.BytesRead += 1
	// read value

	switch aux.Type {
	case 'A':
		value := byte(0)
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Error)

		aux.Value = value
		reader.BytesRead += 1
	case 'c':
		value := int8(0)
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Error)

		aux.Value = value
		aux.Type = 'i'
		reader.BytesRead += 1
	case 'C':
		value := uint8(0)
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Error)

		aux.Value = value
		aux.Type = 'i'
		reader.BytesRead += 1
	case 's':
		value := int16(0)
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Error)

		aux.Value = value
		reader.BytesRead += 2
	case 'S':
		value := uint16(0)
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Error)

		aux.Value = value
		reader.BytesRead += 2
	case 'i':
		value := int32(0)
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Error)

		aux.Value = value
		reader.BytesRead += 4
	case 'I':
		value := uint32(0)
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Error)
		aux.Value = value
		reader.BytesRead += 4
	case 'f':
		value := float32(0)
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Error)

		aux.Value = value
		reader.BytesRead += 4
	case 'd':
		value := float64(0)
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Error)

		aux.Value = value
		reader.BytesRead += 8
	case 'Z':
		var b byte
		var buffer bytes.Buffer
		for {
			reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &b)
			common.ExitIfError(reader.Error)

			reader.BytesRead += 1
			if b == 0 {
				break
			}
			buffer.WriteByte(b)
		}
		aux.Value = buffer.String()
	case 'H':
		var b byte
		var buffer bytes.Buffer
		for {
			reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &b)
			common.ExitIfError(reader.Error)

			reader.BytesRead += 1
			if b == 0 {
				break
			}
			fmt.Fprintf(&buffer, "%X", b)
		}
		aux.Value = buffer.String()
	case 'B':
		var t byte
		var k int32
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &t)
		common.ExitIfError(reader.Error)

		reader.BytesRead += 1
		reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &k)
		common.ExitIfError(reader.Error)

		reader.BytesRead += 4
		switch t {
		case 'c':
			data := make([]int32, k)
			for i = 0; i < int(k); i++ {
				reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &data[i])
				common.ExitIfError(reader.Error)
				reader.BytesRead += 1
			}
			aux.Value = data
		case 'C':
			data := make([]uint8, k)
			for i = 0; i < int(k); i++ {
				reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &data[i])
				common.ExitIfError(reader.Error)

				reader.BytesRead += 1
			}
			aux.Value = data
		case 's':
			data := make([]int16, k)
			for i = 0; i < int(k); i++ {
				reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &data[i])
				common.ExitIfError(reader.Error)

				reader.BytesRead += 2
			}
			aux.Value = data
		case 'S':
			data := make([]uint16, k)
			for i = 0; i < int(k); i++ {
				reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &data[i])
				common.ExitIfError(reader.Error)

				reader.BytesRead += 2
			}
			aux.Value = data
		case 'i':
			data := make([]int32, k)
			for i = 0; i < int(k); i++ {
				reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &data[i])
				common.ExitIfError(reader.Error)

				reader.BytesRead += 4
			}
			aux.Value = data
		case 'I':
			data := make([]uint32, k)
			for i = 0; i < int(k); i++ {
				reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &data[i])
				common.ExitIfError(reader.Error)

				reader.BytesRead += 4
			}
			aux.Value = data
		case 'f':
			data := make([]float32, k)
			for i = 0; i < int(k); i++ {
				reader.Error = binary.Read(reader.gunzip, binary.LittleEndian, &data[i])
				common.ExitIfError(reader.Error)

				reader.BytesRead += 4
			}
			aux.Value = data
		default:
			reader.Error = fmt.Errorf("invalid auxiliary array value type `%c'", t)
			common.ExitIfError(reader.Error)
		}

	default:
		reader.Error = fmt.Errorf("invalid auxiliary value type `%c'", aux.Type)
		common.ExitIfError(reader.Error)
	}
	return aux
}

type CigarByte struct {
	Op  byte
	Len int
}

type cInfo struct {
	Text string
	Len  int32
	NRef int32
}

//allocates memory for new bam header
func MakeHeader() *BamHeader {
	return &BamHeader{ChromSize: make(map[string]int), decoder: &cInfo{Len: 0, NRef: 0}}
}

func bamCigar(bCig []uint32) string {
	return byteCigarToString(ParseCigar(bCig))
}

func ToHeavyCigar(bCig []uint32) []*cigar.Cigar {
	var ans []*cigar.Cigar
	for _, i := range ParseCigar(bCig) {
		ans = append(ans, &cigar.Cigar{RunLength: int64(i.Len), Op: rune(i.Op)})
	}
	return ans
}

func byteCigarToString(cig []CigarByte) string {
	var ans string = ""
	for i := 0; i < len(cig); i++ {
		ans += fmt.Sprintf("%s%d", string(cig[i].Op), cig[i].Len)
	}
	return ans
}

func ParseCigar(bamCigar []uint32) []CigarByte {
	var ans []CigarByte
	types := []byte{'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}
	for i := 0; i < len(bamCigar); i++ {
		n := bamCigar[i] >> 4
		t := types[bamCigar[i]&0xf]
		ans = append(ans, CigarByte{Op: t, Len: int(n)})
	}
	return ans
}
