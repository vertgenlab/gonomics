package bam

import (
	"encoding/binary"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/common"
	"log"
	"strings"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"bytes"
	"bufio"
	"io"
	"fmt"
	"io/ioutil"
	"os"
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
	Note      []string
}

func BamBlockToSam(header *BamHeader, bam *BamData) {
	//var buffer bytes.Buffer
	var ans *sam.SamAln = &sam.SamAln{
		QName: bam.QName,
		Flag: int64(bam.Flag),
		RName: header.Chroms[bam.RName].Name,
		Pos: int64(bam.Pos),
		MapQ: int64(bam.MapQ),
		RNext: "",
		PNext: int64(bam.NextPos),
		TLen: int64(header.Chroms[bam.RName].Size),
		Seq: dna.StringToBases(bamSequence(bam.Seq)),
		Qual: qualToString(bam.Qual),
		//Extra:
	}
	log.Printf("%v\n", sam.SamAlnToString(ans))


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
	if strings.Contains(string(qual), "\n") {
		fmt.Fprintf(writer, "%c", '*')
	} else {
		for i := 0; i < len(qual); i++ {
			fmt.Fprintf(writer, "%c", qual[i]+33)
		}
	}
	
	writer.Flush()
	return buffer.String()
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
	Chroms []*chromInfo.ChromInfo
}

type BamReader struct {
	File 	*os.File
	gunzip 	*Bgzip
	PipeIo  *BgzipPipe
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
		bamHeader.Chroms = append(bamHeader.Chroms, &chromInfo.ChromInfo{Name: strings.Trim(string(reader.PipeIo.Data), "\n\000"), Size: int64(lengthSeq), Order: int64(len(bamHeader.Chroms))})
		//bamHeader.ChromSize[string(reader.PipeIo.Data)] = int(lengthSeq)
	}
	//log.Printf("%s\n", bamHeader.Txt)
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

func Read(filename string) {
	bamFile := NewBamReader(filename)
	defer bamFile.File.Close()
	h := ReadHeader(bamFile)
	binaryData := make(chan *BamData)
	go bamToChan(bamFile, binaryData)
	for each := range binaryData {
		BamBlockToSam(h, each)
	}
}

func bamToChan(reader *BamReader, binaryData chan <- *BamData) {
	var blockSize int32
  	var flagNc    uint32
  	var stats   uint32
  	var err error
  	//binaryData := make(chan *BamData)
  	for {
  		block := &BamData{}
  		buf := bytes.NewBuffer([]byte{})
  		// read block size
  		if err = binary.Read(reader.gunzip, binary.LittleEndian, &blockSize); err != nil {
  			if err == io.EOF {
  				break
      		}
      		common.ExitIfError(err)
  		}
  		//read block data
  		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.RName); err != nil {
  			common.ExitIfError(err)
  		}
  		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Pos); err != nil {
      		common.ExitIfError(err)
      	}
      	if err = binary.Read(reader.gunzip, binary.LittleEndian, &stats); err != nil {
      		common.ExitIfError(err)
      	}
      	block.Bai      = uint16((stats >> 16) & 0xffff)
    	block.MapQ     = uint8 ((stats >>  8) & 0xff)
    	block.RNLength = uint8 ((stats >>  0) & 0xff)
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &flagNc); err != nil {
    		common.ExitIfError(err)
    	}
    	// get Flag and NCigarOp from FlagNc
    	block.Flag     = uint16(flagNc >> 16)
    	block.NCigarOp = uint16(flagNc & 0xffff)
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.LSeq); err != nil {
    		common.ExitIfError(err)
    	}
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextRefID); err != nil {
    		common.ExitIfError(err)
    	}
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextPos); err != nil {
    		common.ExitIfError(err)
    	}
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.TLength); err != nil {
    	 	common.ExitIfError(err)
    	}
    	// parse the read name
    	var b byte
    	for {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &b); err != nil {
    			common.ExitIfError(err)

    		}
    		if b == 0 {
        		block.QName = buf.String()
        		break
      		}
      		buf.WriteByte(b)
    	}
    	var i int
    	// parse cigar block
    	block.Cigar = make(BamCigar, block.NCigarOp)
    	for i := 0; i < int(block.NCigarOp); i++ {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Cigar[i]); err != nil {
    			common.ExitIfError(err)
    		}
    	}
    	 // parse seq
    	block.Seq = make([]byte, (block.LSeq+1)/2)
    	for i = 0; i < int((block.LSeq+1)/2); i++ {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Seq[i]); err != nil {
    			common.ExitIfError(err)
    		}
    	}
    	block.Qual = make([]byte, block.LSeq)
    	for i = 0; i < int(block.LSeq); i++ {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Qual[i]); err != nil {
    			common.ExitIfError(err)
    		}
    	}
    	// read auxiliary data
    	position := 8*4 + int(block.RNLength) + 4*int(block.NCigarOp) + int((block.LSeq + 1)/2) + int(block.LSeq)
    	//for i := 0; position + i < int(blockSize); {

    	//}
    	if _, err = io.CopyN(ioutil.Discard, reader.gunzip, int64(blockSize) - int64(position)); err != nil {
    		common.ExitIfError(err)
    	}
    	//log.Printf("%v\n", block)
    	binaryData <- block
  	}

}

type BamAuxiliary struct {
	Tag   [2]byte
	Value interface{}
}
type BamCigar []uint32

func bamBlocks(bams []*BamData) {
	for i:= 0; i < len(bams);i ++ {
		log.Printf("%v\n", bams[i])
	}
}




type CigarByte struct {
	Op  byte
	Len int
}

type cInfo struct {
	Text    string
	Len 	int32
	NRef 	int32
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


//func bamBlockToSam()
/*
func bamTryAgain(reader *BamReader) {
	h := ReadHeader(reader)
	var blockSize int32
  	var flagNc    uint32
  	var stats   uint32
  	var err error
  	for {
  		block := BamData{}
  		buf := bytes.NewBuffer([]byte{})
  		// read block size
  		if err = binary.Read(reader.gunzip, binary.LittleEndian, &blockSize); err != nil {
  			if err == io.EOF {
  				return
      		}
      		common.ExitIfError(err)
  		}
  		//read block data
  		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.RName); err != nil {
  			common.ExitIfError(err)
  		}
  		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Pos); err != nil {
      		common.ExitIfError(err)
      	}
      	if err = binary.Read(reader.gunzip, binary.LittleEndian, &stats); err != nil {
      		common.ExitIfError(err)
      	}
      	block.Bai      = uint16((stats >> 16) & 0xffff)
    	block.MapQ     = uint8 ((stats >>  8) & 0xff)
    	block.RNLength = uint8 ((stats >>  0) & 0xff)
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &flagNc); err != nil {
    		common.ExitIfError(err)
    	}
    	// get Flag and NCigarOp from FlagNc
    	block.Flag     = uint16(flagNc >> 16)
    	block.NCigarOp = uint16(flagNc & 0xffff)
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.LSeq); err != nil {
    		common.ExitIfError(err)
    	}
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextRefID); err != nil {
    		common.ExitIfError(err)
    	}
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextPos); err != nil {
    		common.ExitIfError(err)
    	}
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.TLength); err != nil {
    	 	common.ExitIfError(err)
    	}
    	// parse the read name
    	var b byte
    	for {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &b); err != nil {
    			common.ExitIfError(err)

    		}
    		if b == 0 {
        		block.QName = buf.String()
        		break
      		}
      		buf.WriteByte(b)
    	}
    	var i int
    	// parse cigar block
    	block.Cigar = make(BamCigar, block.NCigarOp)
    	for i := 0; i < int(block.NCigarOp); i++ {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Cigar[i]); err != nil {
    			common.ExitIfError(err)
    		}
    	}
    	 // parse seq
    	block.Seq = make([]byte, (block.LSeq+1)/2)
    	for i = 0; i < int((block.LSeq+1)/2); i++ {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Seq[i]); err != nil {
    			common.ExitIfError(err)
    		}
    	}
    	block.Qual = make([]byte, block.LSeq)
    	for i = 0; i < int(block.LSeq); i++ {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Qual[i]); err != nil {
    			common.ExitIfError(err)
    		}
    	}
    	// read auxiliary data
    	position := 8*4 + int(block.RNLength) + 4*int(block.NCigarOp) + int((block.LSeq + 1)/2) + int(block.LSeq)
    	//for i := 0; position + i < int(blockSize); {

    	//}
    	if _, err = io.CopyN(ioutil.Discard, reader.gunzip, int64(blockSize) - int64(position)); err != nil {
    		common.ExitIfError(err)
    	}
    	BamBlockToSam(h, &block)
    	//log.Printf("Block Value: %v\n", block)
  	}
}*/
	