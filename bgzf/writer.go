package bgzf

import (
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"encoding/hex"
	"fmt"
	"io"
	"log"
	"math"
)

// magicHexEOF is a particular empty bgzf block that marks the true EOF.
var magicEOF []byte = makeMagicEOF()

// BlockWriter moves data -> compressor -> zipBlock -> file writer
type BlockWriter struct {
	w          io.Writer
	compressor *gzip.Writer
	zipBlock   *Block
}

type Writer struct {
	bw *BlockWriter
	buf bytes.Buffer
}

// NewBlockWriter creates a bgzf writer from any input writer.
func NewBlockWriter(w io.Writer) *BlockWriter {
	var zw BlockWriter
	zw.w = w
	zw.zipBlock = NewBlock()
	zw.compressor = gzip.NewWriter(zw.zipBlock)
	return &zw
}

// NewWriter creates a bgzf writer that writes bytes instead of blocks
func NewWriter(w io.Writer) *Writer {
	var bw Writer
	bw.bw = NewBlockWriter(w)
	return &bw
}

// Write input bytes as a single bgzf block.
func (w *BlockWriter) Write(p []byte) (n int, err error) {
	w.compressor.Reset(w.zipBlock)

	// ******
	// TODO better way around this?
	// The gzip writer tries to write a header smartly at the time of
	// data write, however bgzf requires custom header fields which
	// require knowing the compressed size of the block, which can only
	// be determined after compression. Therefore we cannot use the
	// default gzip header behavior. The code below performs an empty
	// write to make the gzip writer 'think' it has written a header.
	// We then trash the default header and write a custom header below.
	_, err = w.compressor.Write([]byte{}) // trash gzip default header
	if err != nil {
		log.Panic(err)
	}
	w.zipBlock.Reset() // reset written trashed header
	// end of code block for discarding default gzip header
	// ******

	n, err = w.compressor.Write(p) // compress p into write buffer
	if n != len(p) || err != nil {
		fmt.Printf("input %d bytes, wrote %d bytes\n", len(p), n)
		log.Panic(err)
	}
	err = w.compressor.Close()
	if err != nil {
		log.Panic(err)
	}

	w.writeHeader(w.zipBlock.Len())        // write block header to file based on compressed size
	_, err = w.w.Write(w.zipBlock.Bytes()) // write buffer to file
	return
}

// Write writes any number of bytes to a buffer then writes the buffer
// as a block once 64KB of bytes have been stored.
func (w *Writer) Write(p []byte) (n int, err error) {
	w.buf.Write(p)
	n = len(p)
	for w.buf.Len() >= 64000 { // 64KB block size
		_, err = w.bw.Write(w.buf.Next(64000))
	}
	return
}

// Close bgzf writer and add the magic EOF marker to the end of the file.
func (w *BlockWriter) Close() error {
	err := w.compressor.Close()
	if err != nil {
		log.Panic(err)
	}
	_, err = w.w.Write(magicEOF)
	return err
}

// Close bgzf writer and add the magic EOF marker to the end of the file.
func (w *Writer) Close() error {
	var err error
	if w.buf.Len() > 0 {
		_, err = w.bw.Write(w.buf.Bytes())
	}
	closeErr := w.bw.Close()

	if err == nil {
		return closeErr
	}
	return closeErr
}

// writeHeader for custom bgzf header.
func (w *BlockWriter) writeHeader(compSize int) {
	var header [18]byte
	header[0] = 31  // gzip ID 1
	header[1] = 139 // gzip ID 2
	header[2] = 8   // gzip compression method
	header[3] = 4   // gzip flags
	// header[4:8] = mod time uint32
	// header[8] = extra flags (0 for bgzf)
	header[9] = 255                                 // OS = unknown
	binary.LittleEndian.PutUint16(header[10:12], 6) // extra len in bytes
	header[12] = 66                                 // bgzf extra ID 1
	header[13] = 67                                 // bgzf extra ID 2
	binary.LittleEndian.PutUint16(header[14:16], 2) // bgzf extra payload size
	if compSize > math.MaxUint16 {
		log.Panic("block size overflow")
	}
	binary.LittleEndian.PutUint16(header[16:18], uint16(compSize+len(header)-1)) // bgzf compressed block size
	n, err := w.w.Write(header[:])
	if n != 18 || err != nil {
		log.Panic(err)
	}
}

// magicEOF generates the magic EOF byte slice from the hex const.
func makeMagicEOF() []byte {
	// empty gzip block per sam specs
	var magicHexEOF string = "1f8b08040000000000ff0600424302001b0003000000000000000000"
	answer, err := hex.DecodeString(magicHexEOF)
	if err != nil {
		log.Panic(err)
	}
	return answer
}
