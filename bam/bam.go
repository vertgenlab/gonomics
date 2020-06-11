package bam

import (
	"encoding/binary"
	"github.com/vertgenlab/gonomics/common"
	//"github.com/vertgenlab/gonomics/fileio"
	"log"
	//"strings"
	//	"github.com/vertgenlab/gonomics/sam"
	//"github.com/vertgenlab/gonomics/chromInfo"
	"bytes"
	//"bufio"
	"io"
	"io/ioutil"
)

//Note: This is not a Bam struct but palce holders for pointers to decode binary file
//4.2 The Bam format defined in hts specs
type BinAln struct {
	BlockSize int32
	RName     int32
	Pos       int32
	Bai       uint16
	MapQ      uint8
	RNLength  uint8
	Flag      BamFlag
	NCigarOp  uint16
	LSeq      int32
	NextRefID int32
	NextPos   int32
	TLength   int32
	QName     string
	Cigar     []uint32
	Seq       []byte
	Qual      []byte
	Note      []string
}

type BamFlag uint16

func processBamRecord(reader *BamReader) {
	var blockSize int32
	var flagNc uint32
	var binMqNl uint32
	var ans []*BinAln
	//

	reader.PipeIo.Data = make([]byte, blockSize)
	for reader.PipeIo.Debug = binary.Read(bytes.NewBuffer(reader.PipeIo.Data), binary.LittleEndian, &blockSize); reader.PipeIo.Debug == io.EOF && reader.PipeIo.Debug != nil; {

		bai := BinAln{}
		//!ByteErrOF(reader.PipeIo.Debug); {

		buf := bytes.NewBuffer(reader.PipeIo.Data)
		//log.Printf("Bytes read: %d\n",reader.PipeIo.Data)
		//reads the block size
		//reader.PipeIo.Data = make([]byte, blockSize)
		//block data
		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.RName)
		common.ExitIfError(reader.PipeIo.Debug)
		reader.PipeIo.Data = make([]byte, bai.RName)
		//var i, j int = 0, 0
		//	for i = 0; i < int(bai.RName); {
		//		j, reader.PipeIo.Debug = reader.gunzip.Read(reader.PipeIo.Data[i:])
		//		common.ExitIfError(reader.PipeIo.Debug)
		//		i += j
		//	}

		var i int

		//position
		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.Pos)
		common.ExitIfError(reader.PipeIo.Debug)

		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &binMqNl)
		common.ExitIfError(reader.PipeIo.Debug)

		bai.Bai = uint16((binMqNl >> 16) & 0xffff)
		bai.MapQ = uint8((binMqNl >> 8) & 0xff)
		bai.RNLength = uint8((binMqNl >> 0) & 0xff)

		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &flagNc)
		common.ExitIfError(reader.PipeIo.Debug)

		// get Flag and NCigarOp from FlagNc
		bai.Flag = BamFlag(flagNc >> 16)
		bai.NCigarOp = uint16(flagNc & 0xffff)

		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.LSeq)
		common.ExitIfError(reader.PipeIo.Debug)

		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.NextRefID)
		common.ExitIfError(reader.PipeIo.Debug)

		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.NextPos)
		common.ExitIfError(reader.PipeIo.Debug)

		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.TLength)
		common.ExitIfError(reader.PipeIo.Debug)

		var b byte

		for reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &b); reader.PipeIo.Debug != io.EOF && reader.PipeIo.Debug != nil; {

			if b == 0 {
				bai.QName = buf.String()
				break
			}

		}
		bai.Cigar = make([]uint32, bai.NCigarOp)
		for i = 0; i < int(bai.NCigarOp); i++ {
			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.Cigar[i])
			common.ExitIfError(reader.PipeIo.Debug)
		}

		bai.Seq = make([]byte, (bai.LSeq + 1/2))
		for i = 0; i < int((bai.LSeq+1)/2); i++ {
			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.Seq[i])
			common.ExitIfError(reader.PipeIo.Debug)
		}

		bai.Qual = make([]byte, bai.LSeq)
		for i = 0; i < int(bai.LSeq); i++ {
			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.Qual[i])
			common.ExitIfError(reader.PipeIo.Debug)
		}

		position := 8*4 + int(bai.RNLength) + 4*int(bai.NCigarOp) + int((bai.LSeq+1)/2) + int(bai.LSeq)

		for i = 0; position+i < int(blockSize); {
			//reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &auxLen)
			//aux := BamAuxiliary{}

			//reader.gunzip.Read(reader.gunzip, binary.LittleEndian, &bai.TLength)
			//common.ExitIfError(reader.PipeIo.Debug)
			_, reader.PipeIo.Debug = io.CopyN(ioutil.Discard, reader.gunzip, int64(blockSize)-int64(position))
		}
		ans = append(ans, &bai)
	}
	ReadBlock(reader, ans)
}

func bamBlocks(reader *BamReader) {
	//ans := &sam.SamAln{}
	//for r := range reader.Read() {
	//	ans.RName = r.RName
	//	log.Printf("%s\n", ans.Name)
	//}
}

func ReadBlock(reader *BamReader, blocks []*BinAln) {
	//var ans []*sam.SamAln
	var i, j int

	for _, each := range blocks {
		reader.PipeIo.Data = make([]byte, each.BlockSize)
		for i = 0; i < int(each.BlockSize); {
			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &each.BlockSize)
			common.ExitIfError(reader.PipeIo.Debug)
			i += j
		}
		log.Printf("Read %d bytes, value: %s\n", i, reader.PipeIo.Data[:i])
		reader.PipeIo.Data = make([]byte, each.RName)
		for i = 0; i < int(each.RName); {
			j, reader.PipeIo.Debug = reader.gunzip.Read(reader.PipeIo.Data[i:])
			common.ExitIfError(reader.PipeIo.Debug)
			i += j
		}
		log.Printf("Read %d bytes, value: %s\n", i, reader.PipeIo.Data[:i])

	}

}

//type BamFlag uint16

type CigarByte struct {
	Op  byte
	Len int
}

type BamReader struct {
	Bgzip
	PipeIo Pipe
}

type cInfo struct {
	Text string
	Len  int32
	NRef int32
}

func NewReader(filename string) *BamReader {
	var binary *BamReader = &BamReader{}
	var err error
	binary.Reader = NewBgzipReader(filename)

	binary.PipeIo = Pipe{BytesRead: 0, Debug: err}
	return binary
}

type BamHeader struct {
	Txt       string
	ChromSize map[string]int
	decoder   *cInfo
}

//allocates memory for new BamReader header
func MakeHeader() *BamHeader {
	return &BamHeader{ChromSize: make(map[string]int), decoder: &cInfo{Len: 0, NRef: 0}}
}

type Pipe struct {
	Bufio     bytes.Buffer
	Data      []byte
	BytesRead int
	Debug     error
}

func ReadHeader(reader *BamReader) *BamHeader {
	bamHeader := MakeHeader()
	magic := make([]byte, 4)
	magic = BgzipBuffer(&reader.Bgzip, magic)
	if string(magic) != "BamReader\001" {
		log.Fatalf("Not a BamReader file: %s\n", string(reader.PipeIo.Data))
	}
	reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bamHeader.decoder.Len)
	common.ExitIfError(reader.PipeIo.Debug)
	reader.PipeIo.Data = make([]byte, bamHeader.decoder.Len)
	var i, j, k int = 0, 0, 0
	for i = 0; i < int(bamHeader.decoder.Len); {
		j, reader.PipeIo.Debug = reader.gunzip.Read(reader.PipeIo.Data[i:])
		common.ExitIfError(reader.PipeIo.Debug)
		i += j
	}
	bamHeader.Txt = string(reader.PipeIo.Data)

	reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bamHeader.decoder.NRef)
	common.ExitIfError(reader.PipeIo.Debug)

	reader.PipeIo.Data = make([]byte, bamHeader.decoder.NRef)

	var lengthName, lengthSeq int32
	for i = 0; i < int(bamHeader.decoder.NRef); i++ {
		lengthName, lengthSeq = 0, 0
		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &lengthName)
		common.ExitIfError(reader.PipeIo.Debug)

		reader.PipeIo.Data = make([]byte, lengthName)

		for j = 0; j < int(lengthName); {
			k, reader.PipeIo.Debug = reader.gunzip.Read(reader.PipeIo.Data[j:])
			common.ExitIfError(reader.PipeIo.Debug)
			j += k
		}
		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &lengthSeq)
		common.ExitIfError(reader.PipeIo.Debug)

		log.Printf("%s\t%d\n", string(reader.PipeIo.Data), int(lengthSeq))
		bamHeader.ChromSize[string(reader.PipeIo.Data)] = int(lengthSeq)
	}
	return bamHeader
}

func decodeBinaryInt32(reader *BamReader) int {
	var key int32 = 0
	reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &key)
	reader.PipeIo.Data = make([]byte, key)
	var i, j int = 0, 0
	for i = 0; i < int(key); {
		j, reader.PipeIo.Debug = reader.gunzip.Read(reader.PipeIo.Data[i:])
		common.ExitIfError(reader.PipeIo.Debug)
		i += j
	}
	return i
}

func ByteErrOF(err error) bool {
	switch true {
	case err != nil && err != io.EOF:
		common.ExitIfError(err)
		return true
	case err == io.EOF:
		return true
	default:
		return false
	}
}
