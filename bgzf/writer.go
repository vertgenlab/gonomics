package bgzf

import (
	"compress/gzip"
	"encoding/binary"
	"io"
	"log"
	"math"
)

type Writer struct {
	w io.Writer
	compressor *gzip.Writer
	bufCompressor *gzip.Writer
	zipBlock *Block
}

func NewWriter(w io.Writer) Writer {
	var zw Writer
	zw.w = w
	zw.zipBlock = NewBlock()
	zw.compressor = gzip.NewWriter(zw.w)
	zw.bufCompressor = gzip.NewWriter(zw.zipBlock)
	//zw.bufCompressor.Write([]byte("init"))
	//zw.zipBlock.Reset()
	return zw
}

func (w Writer) Write(b []byte) (n int, err error) {
	w.bufCompressor.Reset(w.zipBlock)
	w.zipBlock.Reset()
	_, err = w.bufCompressor.Write(b)
	if err != nil {
		log.Panic(err)
	}
	w.bufCompressor.Flush()

	w.compressor.Reset(w.w)
	w.compressor.Extra = append(w.compressor.Extra, getExtraTag(w.zipBlock.Len())...)
	n, err = w.compressor.Write(b)
	w.compressor.Close()
	return
}

func getExtraTag(len int) []byte {
	if len > math.MaxUint16 {
		log.Panic("buffer is too big")
	}
	answer := make([]byte, 6)
	answer[0] = 'B' // ID[0]
	answer[1] = 'C' // ID[1]
	answer[2] = 2 // payload size uint16[0]
	answer[3] = 0 // payload size uint16[1]
	binary.LittleEndian.PutUint16(answer[4:6], uint16(len) - 1)
	return answer
}

func (w Writer) Close() error {
	w.Write([]byte{})
	//w.compressor.Reset(w.w)
	//w.compressor.Extra = append(w.compressor.Extra, getExtraTag(28)...)
	//w.compressor.Write([]byte{})
	return w.compressor.Close()
}
