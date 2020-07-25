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
	//"github.com/vertgenlab/gonomics/giraf"
	"io"
	//"io/ioutil"
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

	bamR.Stream = &BgZipBuffer{BytesRead: 0, Debug: err}
	return bamR
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
	Stream *BgZipBuffer
}

func ReadHeader(reader *BamReader) *BamHeader {
	bamHeader := MakeHeader()
	magic := make([]byte, 4)
	magic = BgzipBuffer(reader.gunzip, magic)
	if string(magic) != "BAM\001" {
		log.Fatalf("Not a BAM file: %s\n", string(reader.Stream.Data))
	}
	reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bamHeader.decoder.Len)
	common.ExitIfError(reader.Stream.Debug)
	reader.Stream.Data = make([]byte, bamHeader.decoder.Len)
	var i, j, k int = 0, 0, 0
	for i = 0; i < int(bamHeader.decoder.Len); {
		j, reader.Stream.Debug = reader.gunzip.Read(reader.Stream.Data[i:])
		common.ExitIfError(reader.Stream.Debug)
		i += j
	}
	bamHeader.Txt = string(reader.Stream.Data)
	reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bamHeader.decoder.NRef)
	common.ExitIfError(reader.Stream.Debug)
	reader.Stream.Data = make([]byte, bamHeader.decoder.NRef)
	var lengthName, lengthSeq int32
	for i = 0; i < int(bamHeader.decoder.NRef); i++ {
		lengthName, lengthSeq = 0, 0
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &lengthName)
		common.ExitIfError(reader.Stream.Debug)
		reader.Stream.Data = make([]byte, lengthName)
		for j = 0; j < int(lengthName); {
			k, reader.Stream.Debug = reader.gunzip.Read(reader.Stream.Data[j:])
			common.ExitIfError(reader.Stream.Debug)
			j += k
		}
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &lengthSeq)
		common.ExitIfError(reader.Stream.Debug)
		bamHeader.Chroms = append(bamHeader.Chroms, &chromInfo.ChromInfo{Name: strings.Trim(string(reader.Stream.Data), "\n\000"), Size: int64(lengthSeq), Order: int64(len(bamHeader.Chroms))})
		//bamHeader.ChromSize[string(reader.Stream.Data)] = int(lengthSeq)
	}
	//log.Printf("%s\n", bamHeader.Txt)
	return bamHeader
}

func decodeBinaryInt32(reader *BamReader) int {
	var key int32 = 0
	reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &key)
	reader.Stream.Data = make([]byte, key)
	var i, j int = 0, 0
	for i = 0; i < int(key); {
		j, reader.Stream.Debug = reader.gunzip.Read(reader.Stream.Data[i:])
		common.ExitIfError(reader.Stream.Debug)
		i += j
	}
	return i
}

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
	var flagNc uint32
	var stats uint32
	var i, j int
	var b byte
	var block *BamData

	//binaryData := make(chan *BamData)
	for {
		block = &BamData{}
		buf := bytes.NewBuffer([]byte{})
		// read block size
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &blockSize)
		if reader.Stream.Debug == io.EOF {
			close(binaryData)
			break
		}
		common.ExitIfError(reader.Stream.Debug)
		//read block data
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &block.RName)
		common.ExitIfError(reader.Stream.Debug)

		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &block.Pos)
		common.ExitIfError(reader.Stream.Debug)

		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &stats)
		common.ExitIfError(reader.Stream.Debug)

		block.Bai = uint16((stats >> 16) & 0xffff)
		block.MapQ = uint8((stats >> 8) & 0xff)
		block.RNLength = uint8((stats >> 0) & 0xff)

		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &flagNc)
		common.ExitIfError(reader.Stream.Debug)

		// get Flag and NCigarOp from FlagNc
		block.Flag = uint16(flagNc >> 16)
		block.NCigarOp = uint16(flagNc & 0xffff)

		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &block.LSeq)
		common.ExitIfError(reader.Stream.Debug)

		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextRefID)
		common.ExitIfError(reader.Stream.Debug)

		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextPos)
		common.ExitIfError(reader.Stream.Debug)

		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &block.TLength)
		common.ExitIfError(reader.Stream.Debug)

		// parse the read name

		for {
			reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &b)
			common.ExitIfError(reader.Stream.Debug)

			if b == 0 {
				block.QName = buf.String()
				break
			}
			buf.WriteByte(b)
		}
		// parse cigar block
		block.Cigar = make([]uint32, block.NCigarOp)
		for i = 0; i < int(block.NCigarOp) && reader.Stream.Debug == nil; i++ {
			reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &block.Cigar[i])
			common.ExitIfError(reader.Stream.Debug)

		}
		// parse seq
		block.Seq = make([]byte, (block.LSeq+1)/2)
		for i = 0; i < int((block.LSeq+1)/2); i++ {
			reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &block.Seq[i])
			common.ExitIfError(reader.Stream.Debug)

		}
		block.Qual = make([]byte, block.LSeq)
		for i = 0; i < int(block.LSeq); i++ {
			reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &block.Qual[i])
			common.ExitIfError(reader.Stream.Debug)
		}
		// read auxiliary data
		j = 8*4 + int(block.RNLength) + 4*int(block.NCigarOp) + int((block.LSeq+1)/2) + int(block.LSeq)
		for i = 0; j+i < int(blockSize); {
			block.Aux = append(block.Aux, ReadBamAuxiliary(reader))
			i += reader.Stream.BytesRead
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

type BamAuxiliary struct {
	Tag   [2]byte
	Type  byte
	Value interface{}
}

func ReadBamAuxiliary(reader *BamReader) *BamAuxiliary {
	aux := &BamAuxiliary{}
	var i int
	//var valueType byte = 0
	// number of read bytes
	reader.Stream.BytesRead = 0
	// read data
	reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &aux.Tag[0])
	common.ExitIfError(reader.Stream.Debug)

	reader.Stream.BytesRead += 1
	reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &aux.Tag[1])
	common.ExitIfError(reader.Stream.Debug)

	reader.Stream.BytesRead += 1
	reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &aux.Type)
	common.ExitIfError(reader.Stream.Debug)

	// three bytes read so far
	reader.Stream.BytesRead += 1
	// read value

	switch aux.Type {
	case 'A':
		value := byte(0)
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Stream.Debug)

		aux.Value = value
		reader.Stream.BytesRead += 1
	case 'c':
		value := int8(0)
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Stream.Debug)

		aux.Value = value
		aux.Type = 'i'
		reader.Stream.BytesRead += 1
	case 'C':
		value := uint8(0)
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Stream.Debug)

		aux.Value = value
		aux.Type = 'i'
		reader.Stream.BytesRead += 1
	case 's':
		value := int16(0)
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Stream.Debug)

		aux.Value = value
		reader.Stream.BytesRead += 2
	case 'S':
		value := uint16(0)
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Stream.Debug)

		aux.Value = value
		reader.Stream.BytesRead += 2
	case 'i':
		value := int32(0)
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Stream.Debug)

		aux.Value = value
		reader.Stream.BytesRead += 4
	case 'I':
		value := uint32(0)
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Stream.Debug)
		aux.Value = value
		reader.Stream.BytesRead += 4
	case 'f':
		value := float32(0)
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Stream.Debug)

		aux.Value = value
		reader.Stream.BytesRead += 4
	case 'd':
		value := float64(0)
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.Stream.Debug)

		aux.Value = value
		reader.Stream.BytesRead += 8
	case 'Z':
		var b byte
		var buffer bytes.Buffer
		for {
			reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &b)
			common.ExitIfError(reader.Stream.Debug)

			reader.Stream.BytesRead += 1
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
			reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &b)
			common.ExitIfError(reader.Stream.Debug)

			reader.Stream.BytesRead += 1
			if b == 0 {
				break
			}
			fmt.Fprintf(&buffer, "%X", b)
		}
		aux.Value = buffer.String()
	case 'B':
		var t byte
		var k int32
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &t)
		common.ExitIfError(reader.Stream.Debug)

		reader.Stream.BytesRead += 1
		reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &k)
		common.ExitIfError(reader.Stream.Debug)

		reader.Stream.BytesRead += 4
		switch t {
		case 'c':
			tmp := make([]int32, k)
			for i = 0; i < int(k); i++ {
				reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &tmp[i])
				common.ExitIfError(reader.Stream.Debug)
				reader.Stream.BytesRead += 1
			}
			aux.Value = tmp
		case 'C':
			tmp := make([]uint8, k)
			for i = 0; i < int(k); i++ {
				reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &tmp[i])
				common.ExitIfError(reader.Stream.Debug)

				reader.Stream.BytesRead += 1
			}
			aux.Value = tmp
		case 's':
			tmp := make([]int16, k)
			for i = 0; i < int(k); i++ {
				reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &tmp[i])
				common.ExitIfError(reader.Stream.Debug)

				reader.Stream.BytesRead += 2
			}
			aux.Value = tmp
		case 'S':
			tmp := make([]uint16, k)
			for i = 0; i < int(k); i++ {
				reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &tmp[i])
				common.ExitIfError(reader.Stream.Debug)

				reader.Stream.BytesRead += 2
			}
			aux.Value = tmp
		case 'i':
			tmp := make([]int32, k)
			for i = 0; i < int(k); i++ {
				reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &tmp[i])
				common.ExitIfError(reader.Stream.Debug)

				reader.Stream.BytesRead += 4
			}
			aux.Value = tmp
		case 'I':
			tmp := make([]uint32, k)
			for i = 0; i < int(k); i++ {
				reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &tmp[i])
				common.ExitIfError(reader.Stream.Debug)

				reader.Stream.BytesRead += 4
			}
			aux.Value = tmp
		case 'f':
			tmp := make([]float32, k)
			for i = 0; i < int(k); i++ {
				reader.Stream.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &tmp[i])
				common.ExitIfError(reader.Stream.Debug)

				reader.Stream.BytesRead += 4
			}
			aux.Value = tmp
		default:
			reader.Stream.Debug = fmt.Errorf("invalid auxiliary array value type `%c'", t)
			common.ExitIfError(reader.Stream.Debug)
		}

	default:
		reader.Stream.Debug = fmt.Errorf("invalid auxiliary value type `%c'", aux.Type)
		common.ExitIfError(reader.Stream.Debug)
	}
	return aux
}

//type BamCigar []uint32

func bamBlocks(bams []*BamData) {
	for i := 0; i < len(bams); i++ {
		log.Printf("%v\n", bams[i])
	}
}

type CigarByte struct {
	//cigar.Cigar
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

type BgZipBuffer struct {
	Bufio     bytes.Buffer
	Data      []byte
	BytesRead int
	Debug     error
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
