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

// ByteReader wraps bufio.Reader for efficient byte-level file parsing.
type ByteReader struct {
	*bufio.Reader
	File         *os.File
	internalGzip *pgzip.Reader
	line         []byte
	Buffer       *bytes.Buffer
}

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
	file := MustOpen(filename)

	reader := &ByteReader{
		File:   file,
		Buffer: &bytes.Buffer{},
	}

	if strings.HasSuffix(filename, ".gz") {
		gzReader, err := pgzip.NewReader(file)
		exception.PanicOnErr(err)

		reader.internalGzip = gzReader
		reader.Reader = bufio.NewReader(gzReader)
	} else {
		reader.Reader = bufio.NewReader(file)
	}
	return reader
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
	bw.buf = bw.getBuffer()
	return bw
}

// NewWriter returns a new Writer whose buffer has the default size.
func NewByteWriter(filename string) *ByteWriter {
	file := MustCreate(filename)
	writer := NewByteWriterSize(file, defaultBufSize)

	if strings.HasSuffix(filename, ".gz") {
		writer.internalGzip = pgzip.NewWriter(file)
		writer.Writer = NewByteWriterSize(writer.internalGzip, defaultBufSize)
	}
	// TODO: pgzip Writer to write gzip files
	return writer
}

// ReadLine reads a line into Buffer, indicating if more lines are available.
func ReadLine(reader *ByteReader) (*bytes.Buffer, bool) {
	reader.Buffer.Reset() // Reset buffer for new line reading
	var err error
	for reader.line, err = reader.ReadSlice('\n'); err == nil || err == bufio.ErrBufferFull || err == io.EOF; reader.line, err = reader.ReadSlice('\n') {
		if err != nil && err != bufio.ErrBufferFull {
			if err == io.EOF {
				return reader.Buffer, true // End of file reached
			}
			exception.PanicOnErr(err) // Handle unexpected errors
		}
		// Write the read part to the buffer, handling both partial and complete lines
		if len(reader.line) > 0 && reader.line[len(reader.line)-1] == '\n' {
			// Check for carriage return before newline and adjust accordingly
			if len(reader.line) > 1 && reader.line[len(reader.line)-2] == '\r' {
				_, err = reader.Buffer.Write(reader.line[:len(reader.line)-2])
			} else {
				_, err = reader.Buffer.Write(reader.line[:len(reader.line)-1])
			}
			exception.PanicOnErr(err)
			if err == nil {
				return reader.Buffer, false
			}
		} else {
			_, err = reader.Buffer.Write(reader.line)
			exception.PanicOnErr(err)
		}
	}
	return reader.Buffer, false
}

// CatchErrThrowEOF handles EOF errors silently and panics on others.
func CatchErrThrowEOF(err error) {
	if err == io.EOF {
		return
	} else {
		exception.PanicOnErr(err)
	}
}

// bytesToBuffer returns a pointer to Buffer containing reader.line.
func bytesToBuffer(reader *ByteReader) *bytes.Buffer {
	var err error
	if reader.line[len(reader.line)-2] == '\r' {
		_, err = reader.Buffer.Write(reader.line[:len(reader.line)-2])
	} else {
		_, err = reader.Buffer.Write(reader.line[:len(reader.line)-1])
	}
	exception.PanicOnErr(err)
	return reader.Buffer
}

// StringToIntSlice converts comma-separated strings to int slices.
func StringToIntSlice(line string) []int {
	work := strings.Split(line, ",")
	var sliceSize int = len(work)
	if line[len(line)-1] == ',' {
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
