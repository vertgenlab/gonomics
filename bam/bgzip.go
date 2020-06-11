package bam

import (
	//"log"
	"fmt"
	"io"
	//	"log"
	"bufio"
	"compress/gzip"
	"encoding/binary"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
)

//All this information can be found in the sam format manual and specs
const (
	bgzfIndex = "BC\x02\x00\x00\x00"
	minFrame  = 20 + len(bgzfIndex) // Minimum bgzf header+footer length.
	// Magic EOF block is defined in the SAM specification.
	magicBlock = "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
)

const (
	BlockSize    = 0x0ff00 // The maximum size of an uncompressed input data block.
	MaxBlockSize = 0x10000 // The maximum size of a compressed output block.
)

type BgzipHeader struct {
	SI1   uint8
	SI2   uint8
	SLen  uint16
	BSize uint16
}

type Bgzip struct {
	io.Reader
	gunzip *gzip.Reader
	seeker *io.ReadSeeker
	File   *os.File
	err    error
}

func NewBgzipReader(filename string) *Bgzip {
	var ans *Bgzip = &Bgzip{}
	ans.File = fileio.MustOpen(filename)
	ans.gunzip, ans.err = gzip.NewReader(ans.File)
	ans.Reader = bufio.NewReader(ans.gunzip)
	common.ExitIfError(ans.err)
	return ans
}

func BgzipBuffer(file *Bgzip, data []byte) []byte {
	_, err := file.Reader.Read(data)
	common.ExitIfError(err)
	return data
}

func (reader *Bgzip) GetExtra() (*BgzipHeader, error) {
	if len(reader.gunzip.Header.Extra) != 6 {
		return nil, fmt.Errorf("no extra information available")
	}
	extra := BgzipHeader{}
	extra.SI1 = reader.gunzip.Header.Extra[0]
	extra.SI2 = reader.gunzip.Header.Extra[1]
	extra.SLen = binary.LittleEndian.Uint16(reader.gunzip.Header.Extra[2:4])
	extra.BSize = binary.LittleEndian.Uint16(reader.gunzip.Header.Extra[4:6])
	return &extra, nil
}

func (gz *Bgzip) Seek(offset int64, whence int) (int64, error) {
	var x int64

	x, gz.err = gz.Seek(offset, whence)

	if gz.err == nil {
		whenceStr := "invalid"
		switch whence {
		case io.SeekStart:
			whenceStr = "SeekStart"
		case io.SeekCurrent:
			whenceStr = "SeekCurrent"
		case io.SeekEnd:
			whenceStr = "SeekEnd"
		}
		fmt.Printf("Seek(%d, %s) CurPos:%d\n", offset, whenceStr, x)
	}

	return x, gz.err
}
