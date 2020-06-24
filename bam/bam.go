package bam

import (
	"encoding/binary"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/common"
	"log"
	//"strings"
	//"github.com/vertgenlab/gonomics/sam"
	//"github.com/vertgenlab/gonomics/chromInfo"
	"bytes"
	//"bufio"
	"io"
	"io/ioutil"

	"os"
)
//Note: This is not a bam struct but palce holders for pointers to decode binary file
//4.2 The BAM format defined in hts specs
type BinAln struct {
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


func processBamRecord(reader *BamReader) []*BinAln {
	var blockSize uint32
  	var flagNc    uint32
  	var binMqNl   uint32
  	var ans []*BinAln
  	for {
  		
  		//buf := bytes.NewBuffer(reader.PipeIo.Data)
  		bai := BinAln{}
  		//reads the block size
  		reader.PipeIo.Data = make([]byte, blockSize)
  		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &blockSize)
  		
  		common.ExitIfError(reader.PipeIo.Debug)
  		reader.PipeIo.Data = make([]byte, blockSize)
  		reader.PipeIo.BytesRead = decodeBinaryInt32(reader)
  		

  		//block data
  		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.RName)
  		common.ExitIfError(reader.PipeIo.Debug)
  			
  		//position
  		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.Pos)
  		common.ExitIfError(reader.PipeIo.Debug)
  			
  		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &binMqNl)
  		common.ExitIfError(reader.PipeIo.Debug)

  		bai.Bai      = uint16((binMqNl >> 16) & 0xffff)
    	bai.MapQ     = uint8 ((binMqNl >>  8) & 0xff)
    	bai.RNLength = uint8 ((binMqNl >>  0) & 0xff)

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

  		
  	
  		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &b)
  		common.ExitIfError(reader.PipeIo.Debug)


  		bai.Cigar = make([]uint32, bai.NCigarOp)
  		for i := 0; i < int(bai.NCigarOp); i++ {
  			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.Cigar[i])
  			common.ExitIfError(reader.PipeIo.Debug)

  		}

  		bai.Seq = make([]byte, (bai.LSeq+1/2))
  		for i := 0; i < int((bai.LSeq+1)/2); i++ {
  			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.Seq[i])
  			common.ExitIfError(reader.PipeIo.Debug)
  		}

  		bai.Qual = make([]byte, bai.LSeq)
  		for i := 0; i < int(bai.LSeq); i++ {
  			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.Qual[i])
  			common.ExitIfError(reader.PipeIo.Debug)
  		}

  		position := 8*4 + int(bai.RNLength) + 4*int(bai.NCigarOp) + int((bai.LSeq + 1)/2) + int(bai.LSeq)
  		for i := 0; position + i < int(blockSize); {
  			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.TLength)
  			common.ExitIfError(reader.PipeIo.Debug)
  			_, reader.PipeIo.Debug = io.CopyN(ioutil.Discard, reader.gunzip, int64(blockSize) - int64(position))
  		}
  		ans = append(ans, &bai)
  		log.Printf("%v\n", bai)
  	}
  	
  	//return ans
  	//return nil
}

func bamBlocks(bams []*BinAln) {
	for i:= 0; i < len(bams);i ++ {
		log.Printf("%v\n", bams[i])
	}
}

type BamFlag uint16


type CigarByte struct {
	Op  byte
	Len int
}

type BamReader struct {
	File 	*os.File
	gunzip 	*Bgzip
	PipeIo  *BgzipPipe
}

type cInfo struct {
	Text    string
	Len 	int32
	NRef 	int32
}


func NewBamReader(filename string) *BamReader {
	var bamR *BamReader = &BamReader{}
	var err error
	bamR.File = fileio.MustOpen(filename)
	bamR.gunzip = NewBgzipReader(bamR.File)

	bamR.PipeIo = &BgzipPipe{BytesRead: 0, Debug: err}
	return bamR
}


type BamHeader struct {
	Txt         string
	ChromSize   map[string]int
	decoder     *cInfo
}
//allocates memory for new bam header
func MakeHeader() *BamHeader {
	return &BamHeader{ChromSize: make(map[string]int), decoder: &cInfo{Len: 0, NRef: 0}}
}

type BgzipPipe struct {
	Bufio bytes.Buffer
	Data []byte
	BytesRead int
	Debug error
}

func ReadHeader(reader *BamReader) *BamHeader {

	bamHeader := MakeHeader()
	magic := make([]byte, 4)
	magic = BgzipBuffer(reader.gunzip, magic)
	if string(magic) != "BAM\001" {
		log.Fatalf("Not a BAM file: %s\n", string(reader.PipeIo.Data))
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
		bamHeader.ChromSize[string(reader.PipeIo.Data)] = int(lengthSeq)
	}
	log.Printf("%s\n", bamHeader.Txt)
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

	