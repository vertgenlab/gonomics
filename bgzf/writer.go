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
}

func NewWriter(w io.Writer) Writer {
	return Writer{
		w: w,
		compressor: gzip.NewWriter(w),
	}
}

func (w Writer) Write(b []byte) (n int, err error) {
	w.compressor.Reset(w.w)
	w.compressor.Extra = append(w.compressor.Extra, getExtraTag(b)...)
	return w.compressor.Write(b)
}

func getExtraTag(b []byte) []byte {
	if len(b) > math.MaxUint16 {
		log.Panic("buffer is too big")
	}
	answer := make([]byte, 6)
	answer[0] = 'B' // ID[0]
	answer[1] = 'C' // ID[1]
	answer[2] = 2 // payload size uint16[0]
	answer[3] = 0 // payload size uint16[1]
	binary.LittleEndian.PutUint16(answer[4:6], uint16(len(b)) - 1)
	return answer
}

func (w Writer) Close() error {
	return w.compressor.Close()
}
