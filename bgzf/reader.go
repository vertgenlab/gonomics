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

// BlockReader is the decompressor for BGZF files which is at it's core
// a gzip reader. The BGZF reader retains a *os.File to use for
// seek functionality as well as a *bufio.Reader to implement the
// io.ByteReader interface, required for reading gzip files as
// independent concatenated blocks, rather than a multistream.
type BlockReader struct {
	file         *os.File
	buf          *bufio.Reader
	decompressor *gzip.Reader
}

type Reader struct {
	br 	*BlockReader
	blk          *Block
	intermediate bytes.Buffer
	eof bool
}

// Next returns the next n bytes from the bgzf block.
// Next handles reading new block if more bytes are
// requested than are in the currently store block.
func (r *Reader) Next(n int) []byte {
	// we have enough bytes in intermediate
	if r.intermediate.Len() >= n {
		return r.intermediate.Next(n)
	}

	// we need to go to block for extra bytes
	if r.blk.Len() >= n-r.intermediate.Len() {
		if r.intermediate.Len() == 0 {
			return r.blk.Next(n)
		} else {
			blkLen := r.intermediate.Len()
			return append(r.intermediate.Next(n), r.blk.Next(n-blkLen)...)
		}
	}

	// not enough bytes in the currently stored block
	_, err := r.intermediate.ReadFrom(r.blk)
	if err != nil {
		log.Panic(err)
	}

	err = r.br.ReadBlock(r.blk)
	if err == io.EOF {
		r.eof = true
	}

	// return with fewest appends possible
	if r.intermediate.Len() != 0 {
		intLen := r.intermediate.Len()
		return append(r.intermediate.Next(n), r.blk.Next(n-intLen)...)
	}

	return r.blk.Next(n)
}

// Read retrieves the next len(b) bytes from the reader
// and automatically handles block reading.
func (r *Reader) Read(b []byte) (n int, err error) {
	retrieved := r.Next(len(b))
	n = len(retrieved)
	if r.eof {
		err = io.EOF
	}
	copy(b, retrieved)
	return
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

// NewBlockReader opens a new BGZF reader for the input file
func NewBlockReader(filename string) *BlockReader {
	var r BlockReader
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

	return &r
}

// NewReader opens a new BGZF reader that reads bytes instead of blocks
func NewReader(filename string) *Reader {
	var r Reader
	r.br = NewBlockReader(filename)
	r.blk = NewBlock()
	return &r
}

// Seek moves the reader to the beginning of the block denoted by the
// input offset. A subsequent call to ReadBlock will read the Seeked block.
// The whence input can be io.SeekStart, io.SeekCurrent, or io.SeekEnd.
func (r *BlockReader) Seek(offset int64, whence int) (ret int64, err error) {
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

// Seek wraps the Seek function for BlockReader and resets the reader buffer.
func (r *Reader) Seek(offset int64, whence int) (ret int64, err error) {
	r.blk.Reset()
	r.intermediate.Reset()
	r.eof = false
	return r.br.Seek(offset, whence)
}

// ReadBlock reads the next block into the input Block.
// Returns io.EOF at the end of the file.
func (r *BlockReader) ReadBlock(b *Block) error {
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
func (r *BlockReader) Close() error {
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

// Close the file and gzip reader.
func (r *Reader) Close() error {
	return r.br.Close()
}
