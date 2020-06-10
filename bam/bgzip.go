package bam

import (
	//"log"
  	"io"
    "fmt"
  //	"log"
 	"github.com/vertgenlab/gonomics/common"
  "github.com/vertgenlab/gonomics/fileio"
 	"compress/gzip"
 	"encoding/binary"
 // "bufio"
 )
//All this information can be found in the sam format manual and specs
const (
  bgzfIndex = "BC\x02\x00\x00\x00"
  minFrame  = 20 + len(bgzfIndex) // Minimum bgzf header+footer length.
  // Magic EOF block.
  // The magic block is defined in the SAM specification.
  magicBlock = "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
)

type BgzipHeader struct {
  SI1   uint8
  SI2   uint8
  SLen  uint16
  BSize uint16
}

type Bgzip struct {
  gzip.Reader
}

func OpenBgzipFile(filename string) *Bgzip {
    file := fileio.MustOpen(filename)
    return NewBgzipReader(file)
}

func NewBgzipReader(r io.Reader) *Bgzip {
  reader, err := gzip.NewReader(r)
  common.ExitIfError(err)
  return &Bgzip{Reader: *reader}
}

func BgzipBuffer(file *Bgzip, data []byte) []byte {
  _, err := file.Reader.Read(data)
  common.ExitIfError(err)
  return data
}

func (reader *Bgzip) GetExtra() (*BgzipHeader, error) {
  if len(reader.Header.Extra) != 6 {
    return nil, fmt.Errorf("no extra information available")
  }
  extra := BgzipHeader{}
  extra.SI1   = reader.Header.Extra[0]
  extra.SI2   = reader.Header.Extra[1]
  extra.SLen  = binary.LittleEndian.Uint16(reader.Header.Extra[2:4])
  extra.BSize = binary.LittleEndian.Uint16(reader.Header.Extra[4:6])
  return &extra, nil
}

