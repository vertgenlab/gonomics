package fileio

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
	"sync"
	"unicode/utf8"

	"github.com/klauspost/pgzip"
	"github.com/vertgenlab/gonomics/exception"
)

const (
	defaultBufSize = 4096
)

// ByteReader is a bufio.Reader with internal buffering management avoiding excessive allocations.
type ByteReader struct {
	*bufio.Reader // Embedding *bufio.Reader
	File          *os.File

	// Buffer Management Fields
	buf     []byte
	Buffer  *bytes.Buffer
	bufPool sync.Pool

	// pgzip compression and writing
	internalGzip *pgzip.Reader
}

// ByteWriter provides buffered buffer pooling for enhanced concurrency support and efficient memory management with optional pgzip compression.
type ByteWriter struct {
	io.Writer // Embedding io.Writer

	// Buffer Management Fields
	buf     []byte
	n       int
	flushAt int
	nFlush  int
	bufPool sync.Pool

	// Error Handling
	err    error
	closed bool

	// Synchronization Fields
	mtx                sync.Mutex
	inChunkedWriteMode bool
	noChunkedWrite     *sync.Cond
	notFlushing        *sync.Cond
	chunkedWriter      *sync.Cond

	// pgzip compression and writing
	internalGzip *pgzip.Writer
}

// NewByteReader initializes a ByteReader for given filename, supporting p/gzip.
func NewByteReader(filename string) *ByteReader {
	br := &ByteReader{
		File:   MustOpen(filename),
		Buffer: &bytes.Buffer{},
		bufPool: sync.Pool{
			New: func() interface{} { return make([]byte, defaultBufSize) },
		},
	}
	if IsGzip(br.File) {
		var err error
		br.internalGzip, err = pgzip.NewReader(br.File)
		exception.PanicOnErr(err)
		br.Reader = bufio.NewReader(br.internalGzip)
	} else {
		br.Reader = bufio.NewReader(br.File)
	}
	return br
}

// NewByteWriterSize returns a new ByteWriter with a buffer of at least the specified size.
func NewByteWriterSize(w io.Writer, size int) *ByteWriter {
	if size <= 0 {
		size = defaultBufSize
	}
	if size < utf8.UTFMax {
		size = utf8.UTFMax
	}
	m := new(sync.Mutex)
	bw := &ByteWriter{
		Writer:         w,
		buf:            make([]byte, size),
		flushAt:        2 * size,
		noChunkedWrite: sync.NewCond(m),
		notFlushing:    sync.NewCond(m),
		chunkedWriter:  sync.NewCond(m),
		bufPool: sync.Pool{
			New: func() interface{} { return make([]byte, size) },
		},
	}
	return bw
}

// NewWriter returns a new Writer whose buffer has the default size.
func NewByteWriter(filename string) *ByteWriter {
	file := MustCreate(filename)
	bw := NewByteWriterSize(file, defaultBufSize)

	if strings.HasSuffix(filename, "gz") {
		bw.internalGzip = pgzip.NewWriter(file)
		bw.Writer = NewByteWriterSize(bw.internalGzip, defaultBufSize)
	}
	return bw
}

// ReadLine reads a buf into Buffer, indicating if more lines are available.
func ReadLine(br *ByteReader) (*bytes.Buffer, bool) {
	br.Buffer.Reset() // Reset buffer for new line reading
	var err error
	for br.buf, err = br.ReadSlice('\n'); err == nil || err == bufio.ErrBufferFull || err == io.EOF; br.buf, err = br.ReadSlice('\n') {
		if err != nil && err != bufio.ErrBufferFull {
			if err == io.EOF {
				return br.Buffer, true // End of file reached
			}
			exception.PanicOnErr(err) // Handle unexpected errors
		}
		// Write the read part to the buffer, handling both partial and complete lines
		if len(br.buf) > 0 && br.buf[len(br.buf)-1] == '\n' {
			// Check for carriage return before newline and adjust accordingly
			if len(br.buf) > 1 && br.buf[len(br.buf)-2] == '\r' {
				_, err = br.Buffer.Write(br.buf[:len(br.buf)-2])
			} else {
				_, err = br.Buffer.Write(br.buf[:len(br.buf)-1])
			}
			exception.PanicOnErr(err)
			if err == nil {
				return br.Buffer, false
			}
		} else {
			_, err = br.Buffer.Write(br.buf)
			exception.PanicOnErr(err)
		}
	}
	return br.Buffer, false
}

// CatchErrThrowEOF handles EOF errors silently and panics on others.
func CatchErrThrowEOF(err error) {
	if err == io.EOF {
		return
	} else {
		exception.PanicOnErr(err)
	}
}

// StringToIntSlice converts comma-separated strings to int slices.
func StringToIntSlice(buf string) []int {
	work := strings.Split(buf, ",")
	var sliceSize int = len(work)
	if buf[len(buf)-1] == ',' {
		sliceSize--
	}
	var answer []int = make([]int, sliceSize)
	var err error
	for i := 0; i < sliceSize; i++ {
		answer[i], err = strconv.Atoi(work[i])
		exception.PanicOnErr(err)
	}
	return answer
}

// IntSliceToString converts int slices to comma-separated strings.
func IntSliceToString(nums []int) string {
	ans := strings.Builder{}
	ans.Grow(2 * len(nums))
	for i := 0; i < len(nums); i++ {
		ans.WriteString(fmt.Sprintf("%d", nums[i]))
		ans.WriteByte(',')
	}
	return ans.String()
}

// IntToString converts an int to a string.
func IntToString(i int) string {
	return fmt.Sprintf("%d", i)
}
