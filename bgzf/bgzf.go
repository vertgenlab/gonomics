package bgzf

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"io"
	"log"
	"math"
	"os"
)

type Reader struct {
	file         *os.File
	buf          *bufio.Reader
	decompressor *gzip.Reader
	offset       int64
}

type Block struct {
	bytes.Buffer // cap == 2^16
}

func NewBlock() *Block {
	b := new(Block)
	b.Grow(math.MaxUint16)
	return b
}

type Offset struct {
	Compressed   int64  // offset in compressed file
	Uncompressed uint16 // offset in uncompressed block
}

func NewReader(filename string) Reader {
	var r Reader
	var err error

	r.file, err = os.Open(filename)
	if err != nil {
		log.Panic(err)
	}

	r.buf = bufio.NewReader(r.file) // for ReadByte method

	r.decompressor, err = gzip.NewReader(r.buf)
	if err != nil {
		log.Panic(err)
	}

	r.decompressor.Multistream(false)

	return r
}

func (r Reader) SeekBlock(off Offset) {
	_, err := r.file.Seek(off.Compressed, io.SeekStart)
	if err != nil {
		log.Panic(err)
	}
	r.buf.Reset(r.file)
	err = r.decompressor.Reset(r.buf)
	if err != nil {
		log.Panic(err)
	}
}

func (r Reader) ReadBlock(b *Block) error {
	r.decompressor.Multistream(false) // parse each bgzf block separately
	b.Reset()                         // remove anything left in the input block
	//fmt.Printf("Extra: %v\n", r.decompressor.Extra) // DEBUG
	_, err := b.ReadFrom(r.decompressor) // read into block from gzip reader
	if err != nil {
		log.Panic(err)
	}

	// reset gzip reader to move to next block
	err = r.decompressor.Reset(r.buf) // only returns io.EOF at true EOF
	return err
}

func (r Reader) Close() error {
	err := r.file.Close()
	if err != nil {
		return err
	}
	err = r.decompressor.Close()
	if err != nil {
		log.Panic(err)
	}
	return err
}
