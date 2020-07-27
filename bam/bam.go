//Package bam is used to process binary alignment files and and decode data into human readable sam text.
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

//BamReader contains data fields used to read and process binary file.
type BamReader struct {
	File      *os.File
	gunzip    *Bgzip
	data      []byte
	bytesRead int
	error     error
}

//binaryDecoder: Uses the BAM format 4.2, defined in hts specs. Note: This is not a bam struct but place holders for pointers to decode binary file
type binaryDecoder struct {
	RefID     int32
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
	Aux       []*BamAux
	BlockSize int32
}

//BamHeader is a data structure containing header information
//parsed from the binary alignment file as well as a decoder
//used to convert the bam to sam
type BamHeader struct {
	Txt       string
	ChromSize map[string]int
	decoder   *chrInfo
	Chroms    []*chromInfo.ChromInfo
}

//chrInfo is used to decode fields contained in the header
type chrInfo struct {
	text   string
	length int32
	numRef int32
}

//BamAux is a struct that organizes the extra tags at the end of sam/bam records
type BamAux struct {
	Tag   [2]byte
	Type  byte
	Value interface{}
}

//CigarByte is a light weight cigar stuct
//that stores the runlength as an int (not int64) and Op as a byte.
type CigarByte struct {
	Len int
	Op  byte
}

//Read will process a bam file and return a slice of sam records that were decoded from binary.
func Read(filename string) (*BamHeader, []*sam.SamAln) {
	bamFile := NewBamReader(filename)
	defer bamFile.File.Close()
	h := ReadHeader(bamFile)
	binaryData := make(chan *binaryDecoder)
	var ans []*sam.SamAln
	go BamToChannel(bamFile, binaryData)
	for each := range binaryData {
		ans = append(ans, BamBlockToSam(h, each))
	}
	return h, ans
}

//NewBamReader is similar to fileio.EasyOpen/fileio.EasyReader
//which will allocate memory for the struct fields
//and is ready to start processing bam lines after calling this function.
func NewBamReader(filename string) *BamReader {
	var bamR *BamReader = &BamReader{}
	var err error
	bamR.File = fileio.MustOpen(filename)
	bamR.gunzip = NewBgzipReader(bamR.File)
	bamR.error = err
	return bamR
}

//ReadHeader will take a BamReader structure as an input
//performs a quick check to make sure the binary file is a valid bam
//then process header lines and returns a BamHeader (similar to samHeader).
func ReadHeader(reader *BamReader) *BamHeader {
	bamHeader := MakeHeader()
	magic := make([]byte, 4)
	magic = BgzipBuffer(reader.gunzip, magic)
	if string(magic) != "BAM\001" {
		log.Fatalf("Not a BAM file: %s\n", string(reader.data))
	}
	reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &bamHeader.decoder.length)
	common.ExitIfError(reader.error)
	reader.data = make([]byte, bamHeader.decoder.length)
	var i, j, k int = 0, 0, 0
	for i = 0; i < int(bamHeader.decoder.length); {
		j, reader.error = reader.gunzip.Read(reader.data[i:])
		common.ExitIfError(reader.error)
		i += j
	}
	bamHeader.Txt = string(reader.data)
	reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &bamHeader.decoder.numRef)
	common.ExitIfError(reader.error)
	reader.data = make([]byte, bamHeader.decoder.numRef)
	var lengthName, lengthSeq int32
	for i = 0; i < int(bamHeader.decoder.numRef); i++ {
		lengthName, lengthSeq = 0, 0
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &lengthName)
		common.ExitIfError(reader.error)
		reader.data = make([]byte, lengthName)
		for j = 0; j < int(lengthName); {
			k, reader.error = reader.gunzip.Read(reader.data[j:])
			common.ExitIfError(reader.error)
			j += k
		}
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &lengthSeq)
		common.ExitIfError(reader.error)
		bamHeader.Chroms = append(bamHeader.Chroms, &chromInfo.ChromInfo{Name: strings.Trim(string(reader.data), "\n\000"), Size: int64(lengthSeq), Order: int64(len(bamHeader.Chroms))})
	}
	return bamHeader
}

//BamToChannel is a goroutine that will use the binary reader to decode bam records
//and send them off to a channel that could be processed into a sam record further downstream.
func BamToChannel(reader *BamReader, binaryData chan<- *binaryDecoder) {
	var blockSize int32
	var bitFlag uint32
	var stats uint32
	var i, j int
	var b byte
	var block *binaryDecoder
	//This loop will attempt to decode a bam block, which is equivalent to one sam line in this case
	for {
		block = &binaryDecoder{}
		buf := bytes.NewBuffer([]byte{})
		// read block size
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &blockSize)
		if reader.error == io.EOF {
			close(binaryData)
			break
		}
		common.ExitIfError(reader.error)
		//read block data
		block.BlockSize = blockSize
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &block.RefID)
		common.ExitIfError(reader.error)

		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &block.Pos)
		common.ExitIfError(reader.error)

		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &stats)
		common.ExitIfError(reader.error)

		block.Bai = uint16((stats >> 16) & 0xffff)
		block.MapQ = uint8((stats >> 8) & 0xff)
		block.RNLength = uint8((stats >> 0) & 0xff)

		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &bitFlag)
		common.ExitIfError(reader.error)

		// get Flag and NCigarOp from bitFlag
		block.Flag = uint16(bitFlag >> 16)
		block.NCigarOp = uint16(bitFlag & 0xffff)

		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &block.LSeq)
		common.ExitIfError(reader.error)

		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextRefID)
		common.ExitIfError(reader.error)

		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextPos)
		common.ExitIfError(reader.error)

		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &block.TLength)
		common.ExitIfError(reader.error)

		// parse the read name

		for {
			reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &b)
			common.ExitIfError(reader.error)

			if b == 0 {
				block.QName = buf.String()
				break
			}
			buf.WriteByte(b)
		}
		// parse cigar block
		block.Cigar = make([]uint32, block.NCigarOp)
		for i = 0; i < int(block.NCigarOp) && reader.error == nil; i++ {
			reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &block.Cigar[i])
			common.ExitIfError(reader.error)

		}
		// parse seq
		block.Seq = make([]byte, (block.LSeq+1)/2)
		for i = 0; i < int((block.LSeq+1)/2); i++ {
			reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &block.Seq[i])
			common.ExitIfError(reader.error)

		}
		block.Qual = make([]byte, block.LSeq)
		for i = 0; i < int(block.LSeq); i++ {
			reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &block.Qual[i])
			common.ExitIfError(reader.error)
		}
		// read auxiliary data
		j = 8*4 + int(block.RNLength) + 4*int(block.NCigarOp) + int((block.LSeq+1)/2) + int(block.LSeq)
		for i = 0; j+i < int(blockSize); {
			block.Aux = append(block.Aux, decodeAuxiliary(reader))
			i += reader.bytesRead
		}
		binaryData <- block
	}
}

//BamBlockToSam is a function that will convert a decoded
//(already processed binary structure) to a human readable sam data.
func BamBlockToSam(header *BamHeader, bam *binaryDecoder) *sam.SamAln {
	return &sam.SamAln{
		QName: bam.QName,
		Flag:  int64(bam.Flag),
		RName: header.Chroms[bam.RefID].Name,
		Pos:   int64(bam.Pos) + 1,
		MapQ:  int64(bam.MapQ),
		Cigar: ToHeavyCigar(bam.Cigar),
		RNext: setRNext(header, bam),
		PNext: int64(bam.NextPos) + 1,
		TLen:  int64(bam.TLength),
		Seq:   dna.StringToBases(BamSeq(bam.Seq)),
		Qual:  qualToString(bam.Qual),
		Extra: auxToString(bam.Aux),
	}
}

//setRNext will process the reference name of the mate pair alignment,
//if the alignment is on the same fragment, then will set to "=".
func setRNext(header *BamHeader, bam *binaryDecoder) string {
	if bam.NextRefID == bam.RefID {
		return "="
	} else if bam.NextRefID > 0 {
		return header.Chroms[bam.NextRefID].Name //fmt.Sprintf("Chrom: %d", bam.NextRefID)
	} else {
		return "*"
	}
}

//BamSeq will convert raw bytes to a string which can be converted to dna.Base.
//TODO: Look into converting bytes straight to dna.Base
func BamSeq(seq []byte) string {
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

//qualToString will convert the aul bytes into a readable string.
func qualToString(qual []byte) string {
	var buffer bytes.Buffer
	writer := bufio.NewWriter(&buffer)
	fmt.Fprintf(writer, "%s", string(formatQual(qual)))
	writer.Flush()
	return buffer.String()
}

//formatQual is a helper function that will add the 33 offset to the ASCII values or set '*' if the qual scores do not exist in the bam.
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

//axtToString will convert the sam/bam auxiliary struct into a human readable string.
func auxToString(aux []*BamAux) string {
	var ans []string
	for i := 0; i < len(aux); i++ {
		ans = append(ans, fmt.Sprintf("%c%c:%c:%v", aux[i].Tag[0], aux[i].Tag[1], aux[i].Type, aux[i].Value))
	}
	return strings.Join(ans, "\t")
}

//decodeAuxiliary will use the bam reader struct to decode binary text to sam auxilary fields.
//In giraf this is what we are calling notes.
//TODO: Look to synchronize auxilary and notes.
func decodeAuxiliary(reader *BamReader) *BamAux {
	aux := &BamAux{}
	var i int
	// number of read bytes
	reader.bytesRead = 0
	// read data
	reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &aux.Tag[0])
	common.ExitIfError(reader.error)

	reader.bytesRead += 1
	reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &aux.Tag[1])
	common.ExitIfError(reader.error)

	reader.bytesRead += 1
	reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &aux.Type)
	common.ExitIfError(reader.error)

	// three bytes read so far
	reader.bytesRead += 1
	// read value

	switch aux.Type {
	case 'A':
		value := byte(0)
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.error)

		aux.Value = value
		reader.bytesRead += 1
	case 'c':
		value := int8(0)
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.error)

		aux.Value = value
		aux.Type = 'i'
		reader.bytesRead += 1
	case 'C':
		value := uint8(0)
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.error)

		aux.Value = value
		aux.Type = 'i'
		reader.bytesRead += 1
	case 's':
		value := int16(0)
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.error)
		aux.Type = 'i'
		aux.Value = value
		reader.bytesRead += 2
	case 'S':
		value := uint16(0)
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.error)
		aux.Type = 'i'
		aux.Value = value
		reader.bytesRead += 2
	case 'i':
		value := int32(0)
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.error)
		aux.Type = 'i'
		aux.Value = value
		reader.bytesRead += 4
	case 'I':
		value := uint32(0)
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.error)
		aux.Value = value
		reader.bytesRead += 4
	case 'f':
		value := float32(0)
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.error)

		aux.Value = value
		reader.bytesRead += 4
	case 'd':
		value := float64(0)
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &value)
		common.ExitIfError(reader.error)

		aux.Value = value
		reader.bytesRead += 8
	case 'Z':
		var b byte
		var buffer bytes.Buffer
		for {
			reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &b)
			common.ExitIfError(reader.error)
			reader.bytesRead += 1
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
			reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &b)
			common.ExitIfError(reader.error)
			reader.bytesRead += 1
			if b == 0 {
				break
			}
			fmt.Fprintf(&buffer, "%X", b)
		}
		aux.Value = buffer.String()
	case 'B':
		var t byte
		var k int32
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &t)
		common.ExitIfError(reader.error)

		reader.bytesRead += 1
		reader.error = binary.Read(reader.gunzip, binary.LittleEndian, &k)
		common.ExitIfError(reader.error)

		reader.bytesRead += 4
		switch t {
		case 'c':
			data := make([]int32, k)
			for i, reader.error = 0, binary.Read(reader.gunzip, binary.LittleEndian, &data[i]); i < int(k) && reader.error == nil; i++ {
				reader.bytesRead += 1
			}
			aux.Value = data
		case 'C':
			data := make([]uint8, k)
			for i, reader.error = 0, binary.Read(reader.gunzip, binary.LittleEndian, &data[i]); i < int(k) && reader.error == nil; i++ {
				reader.bytesRead += 1
			}
			aux.Value = data
		case 's':
			data := make([]int16, k)
			for i, reader.error = 0, binary.Read(reader.gunzip, binary.LittleEndian, &data[i]); i < int(k) && reader.error == nil; i++ {
				reader.bytesRead += 2
			}
			aux.Value = data
		case 'S':
			data := make([]uint16, k)
			for i, reader.error = 0, binary.Read(reader.gunzip, binary.LittleEndian, &data[i]); i < int(k) && reader.error == nil; i++ {
				reader.bytesRead += 2
			}
			aux.Value = data
		case 'i':
			data := make([]int32, k)
			for i, reader.error = 0, binary.Read(reader.gunzip, binary.LittleEndian, &data[i]); i < int(k) && reader.error == nil; i++ {
				reader.bytesRead += 4
			}
			aux.Value = data
		case 'I':
			data := make([]uint32, k)
			for i, reader.error = 0, binary.Read(reader.gunzip, binary.LittleEndian, &data[i]); i < int(k) && reader.error == nil; i++ {
				reader.bytesRead += 4
			}
			aux.Value = data
		case 'f':
			data := make([]float32, k)
			for i, reader.error = 0, binary.Read(reader.gunzip, binary.LittleEndian, &data[i]); i < int(k) && reader.error == nil; i++ {
				reader.bytesRead += 4
			}
			aux.Value = data
		default:
			log.Fatalf("Error: encountered unknown auxiliary array value type %c...\n", t)
		}

	default:
		log.Fatalf("Error: Found invalid auxiliary value type %c...\n", aux.Type)
	}
	return aux
}

//MakeHeader allocates memory for new bam header.
func MakeHeader() *BamHeader {
	return &BamHeader{ChromSize: make(map[string]int), decoder: &chrInfo{length: 0, numRef: 0}}
}

//ToHeavyCigar will convert the cigar byte to the cigar struct we are using in the cigar package.
func ToHeavyCigar(bCig []uint32) []*cigar.Cigar {
	var ans []*cigar.Cigar
	for _, i := range ParseCigar(bCig) {
		ans = append(ans, &cigar.Cigar{RunLength: int64(i.Len), Op: rune(i.Op)})
	}
	return ans
}

//ParseCigar will parse and convert a slice of uint32 and return a slice of cigar bytes.
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

//ByteCigarToString will process the cigar byte struct and parse and/or convert the data into a string.
func ByteCigarToString(cig []CigarByte) string {
	var ans string = ""
	for i := 0; i < len(cig); i++ {
		ans += fmt.Sprintf("%s%d", string(cig[i].Op), cig[i].Len)
	}
	return ans
}

/*
// TODO: Plans for the new SAM/BAM record.
type Record struct {
	QName      string
	Ref       string
	Pos       int
	MapQ      byte
	Cigar     []uint32
	Flags     uint16
	MateRef   string
	MatePos   int
	TempLen   int
	Seq       []dna.Base
	Qual      []byte
	AuxFields AuxFields
}*/
