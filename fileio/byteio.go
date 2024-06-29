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
	"unsafe"

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
	err error

	// Buffer Management Fields
	buf     []byte
	Buffer  *bytes.Buffer
	bufPool sync.Pool

	// pgzip compression and writing
	internalGzip *pgzip.Reader
}

type ByteWriter struct {
	// Protects all internal state
	mtx *sync.Mutex

	err error
	buf []byte
	n   int
	wr  *bufio.Writer
	bufPool sync.Pool

	// Whether a chunked write is in flight and Writer is in "chunked write mode"
	inChunkedWriteMode bool
	// Condition variable for all goroutines waiting on a chunked write
	noChunkedWrite *sync.Cond

	// Fields below are only used by flush() (and marginally by Reset()) to
	// provide correct serialization of writes.

	// Number of bytes currently being flushed from the start of buf. Only
	// non-zero while a flush is in progress. Always <= n.
	nFlush int
	// Condition variable for all goroutines waiting to flush()
	notFlushing *sync.Cond
	// Condition variable that the only chunked writer is waiting on to flush()
	chunkedWriter *sync.Cond

	// Automatically flush when n > flushAt
	flushAt int

	internalGzip *pgzip.Writer
	gzipPool sync.Pool
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
		br.internalGzip, br.err = pgzip.NewReader(br.File)
		exception.PanicOnErr(br.err)
		br.Reader = bufio.NewReader(br.internalGzip)
	} else {
		br.Reader = bufio.NewReader(br.File)
	}
	br.buf = br.getBuffer()
	return br
}

// NewWriter returns a new Writer whose buffer has the default size.
func NewByteWriter(filename string) *ByteWriter {
	file := MustCreate(filename)
	bw := NewByteWriterSize(file, defaultBufSize)

	if strings.HasSuffix(filename, "gz") {
		bw.gzipPool = sync.Pool{ 
			New: func() interface{} { return pgzip.NewWriter(nil) },
		}
		bw.internalGzip = bw.gzipPool.Get().(*pgzip.Writer)
        bw.internalGzip.Reset(file)
	}
	return bw
}


// NewWriterSize returns a new Writer whose buffer has at least the specified
// size. If the argument io.Writer is already a Writer with large enough size,
// it returns the underlying Writer.
func NewByteWriterSize(w io.Writer, size int) *ByteWriter {
	if size <= 0 {
		size = defaultBufSize
	}
	if size < utf8.UTFMax {
		size = utf8.UTFMax
	}
	m := new(sync.Mutex)
	bw:= &ByteWriter{
		mtx:            m,
		bufPool: sync.Pool{
			New: func() interface{} { return make([]byte, defaultBufSize) },
		},
		wr:             bufio.NewWriter(w),
		noChunkedWrite: sync.NewCond(m),
		notFlushing:    sync.NewCond(m),
		chunkedWriter:  sync.NewCond(m),
		flushAt:        2 * size,
	}
	bw.buf = bw.getBuffer()
	return bw
}

// ReadLine reads a buf into Buffer, indicating if more lines are available.
func ReadLine(br *ByteReader) (*bytes.Buffer, bool) {
	br.Buffer.Reset() // Reset buffer for new line reading
	var err error
	if br.buf != nil {
		br.putBuffer(br.buf)
	}
	br.buf = br.getBuffer()
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

func NextByteLine(br *ByteReader) ([]byte, bool) {
    line := br.getBuffer()
    defer br.putBuffer(line)

    for br.buf, br.err = br.ReadSlice('\n'); br.err == nil || br.err == bufio.ErrBufferFull || br.err == io.EOF; br.buf, br.err = br.ReadSlice('\n')  {
        if br.err != nil && br.err != bufio.ErrBufferFull {
            if br.err == io.EOF {
                if len(br.buf) > 0 {
                    line = append(line, br.buf...) // Append last partial line
                }
                return line, true // EOF reached (or no more complete lines)
            } else {
                exception.PanicOnErr(br.err)
                return nil, true // Signal end of reading due to error
            }
        }
    }
	return bytes.TrimRight(append(line, br.buf...), "\r\n"), false
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

func StringToBytes(s string) []byte {
	p := unsafe.StringData(s)
	b := unsafe.Slice(p, len(s))
	return b
}