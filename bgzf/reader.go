package bgzf

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"log"
	"math"
	"os"
)

// Reader is the decompressor for BGZF files which is at it's core
// a gzip reader. The BGZF reader retains a *os.File to use for
// seek functionality as well as a *bufio.Reader to implement the
// io.ByteReader interface, required for reading gzip files as
// independent concatenated blocks, rather than a multistream.
type Reader struct {
	file         *os.File
	buf          *bufio.Reader
	decompressor *gzip.Reader
}

// Block represents a single BGZF block stores as a bytes.Buffer
// initialized to have a capacity of 65536 (max uint16).
type Block struct {
	bytes.Buffer // cap == 2^16
}

// NewBlock returns a pointer to a bytes buffer with capacity 65536.
func NewBlock() *Block {
	b := new(Block)
	b.Grow(math.MaxUint16)
	return b
}

// Offset for a location in the BGZF file consisting of the compressed
// offset (the number of bytes before the beginning of the desired block
// starting from the beginning of the file), and and the uncompressed
// offset (the number of bytes from the desired record starting from the
// beginning of the uncompressed block)
type Offset struct {
	Compressed   int64  // offset in compressed file
	Uncompressed uint16 // offset in uncompressed block
}

// NewReader opens a new BGZF reader for the input file
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

	return r
}

// Seek moves the reader to the beginning of the block denoted by the
// input offset. A subsequent call to ReadBlock will read the Seeked block.
// The whence input can be io.SeekStart, io.SeekCurrent, or io.SeekEnd.
func (r Reader) Seek(offset int64, whence int) (ret int64, err error) {
	ret, err = r.file.Seek(offset, whence)
	if err != nil {
		log.Panic(err)
	}
	r.buf.Reset(r.file)
	resetErr := r.decompressor.Reset(r.buf)
	if resetErr != nil {
		log.Panic(resetErr)
	}
	return
}

// ReadBlock reads the next block into the input Block.
// Returns io.EOF at the end of the file.
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

// Close the file and gzip reader.
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
