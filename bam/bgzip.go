package bam

import (
	"compress/gzip"
	"encoding/binary"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
)

//Bgzip is a struct that contains an embedded gzip.Reader
//TODO: Seek functionality is not implemented and if we move towards that direction, bgzip could become its own package.
type Bgzip struct {
	gzip.Reader
}

//BgzipHeader contains the extra information that is at the top of the bgzip file that is used to seek to a specific line in the bgzip file.
type BgzipHeader struct {
	SI1   uint8
	SI2   uint8
	SLen  uint16
	BSize uint16
}

//OpenBgzipFile is a wrapper that will open a file and return a gzip reader that can further process data in the file.
func OpenBgzipFile(filename string) *Bgzip {
	file := fileio.MustOpen(filename)
	return NewBgzipReader(file)
}

//NewBgzipReader will take an io.Reader and return a bgzip struct that contains an embedded gzip.Reader
func NewBgzipReader(r io.Reader) *Bgzip {
	reader, err := gzip.NewReader(r)
	common.ExitIfError(err)
	return &Bgzip{Reader: *reader}
}

//BgzipBuffer will gunzip text, catch any possible errors, and return a slice of bytes.
func BgzipBuffer(file *Bgzip, data []byte) []byte {
	_, err := file.Reader.Read(data)
	common.ExitIfError(err)
	return data
}

//GetExtra is a method implemented to retrieve header information used to seek lines in a file.
func (reader *Bgzip) GetExtra() (*BgzipHeader, error) {
	if len(reader.Header.Extra) != 6 {
		return nil, fmt.Errorf("no extra information available")
	}
	extra := BgzipHeader{}
	extra.SI1 = reader.Header.Extra[0]
	extra.SI2 = reader.Header.Extra[1]
	extra.SLen = binary.LittleEndian.Uint16(reader.Header.Extra[2:4])
	extra.BSize = binary.LittleEndian.Uint16(reader.Header.Extra[4:6])
	return &extra, nil
}
