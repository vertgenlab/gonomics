package fileio

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"github.com/vertgenlab/gonomics/common"
	"io"
	"log"
	"os"
	"strings"
	"sync"
)

const (
	defaultBufSize = 4096
)

// SimpleReader implements the io.Reader interface by providing
// the Read(b []byte) method. The struct contains an embedded *bufio.Reader
// and a pointer to os.File for closeure when reading is complete.
type SimpleReader struct {
	*bufio.Reader
	File   *os.File
	Pool   *sync.Pool
	Buffer *bytes.Buffer
}

type Line struct {
	Pool []byte
}

func NewLine() *Line {
	return &Line{
		Pool: make([]byte, defaultBufSize),
	}
}

// Read reads data into p. It returns the number of bytes read into p.
func (reader *SimpleReader) Read(b []byte) (n int, err error) {
	return reader.Read(b)
}

// NewSimpleReader will process a given file and performs error handling if an error occurs.
// SimpleReader will prcoess gzipped files accordinging by performing a check on the suffix
// of the provided file.
func NewSimpleReader(filename string) *SimpleReader {
	var pool = sync.Pool{
		// New creates an object when the pool has nothing available to return.
		// New must return an interface{} to make it flexible. You have to cast
		// your type after getting it.
		New: func() interface{} {
			// Pools often contain things like *bytes.Buffer, which are
			// temporary and re-usable. In this case we have a pointer to a slice of bytes.
			return NewLine()
		},
	}
	var answer SimpleReader = SimpleReader{
		File:   MustOpen(filename),
		Pool:   &pool,
		Buffer: &bytes.Buffer{},
	}
	switch true {
	case strings.HasSuffix(filename, ".gz"):
		gzipReader, err := gzip.NewReader(answer.File)
		common.ExitIfError(err)
		answer.Reader = bufio.NewReader(gzipReader)
	default:
		answer.Reader = bufio.NewReader(answer.File)
	}
	return &answer
}

// ReadLine will return a bytes.Buffer pointing to the internal slice of bytes. Provided this function is called within a loop,
// the function will read one line at a time, and return bool to continue reading. Important to note the buffer return points to
// the internal slice belonging to the reader, meaning the slice will be overridden if the data is not copied. Please be aware the
// reader will call close on the file once the reader encounters EOF.
func ReadLine(reader *SimpleReader) (*bytes.Buffer, bool) {
	var err error
	curr := reader.Pool.Get().(*Line)
	defer reader.Pool.Put(curr)
	curr.Pool = curr.Pool[:0]
	curr.Pool, err = reader.ReadSlice('\n')
	if err == nil {
		if curr.Pool[len(curr.Pool)-1] == '\n' {
			reader.Buffer.Reset()
			reader.Buffer.Write(curr.Pool[:len(curr.Pool)-1])
			return reader.Buffer, false
		} else {
			log.Fatalf("Error: end of line did not end with a white space character...\n")
		}
	}
	CatchErrThrowEOF(err)
	reader.Close()
	return nil, true
}

// CatchErrThrowEOF will silently handles and throws the EOF error and will log and exit any other errors.
func CatchErrThrowEOF(err error) {
	if err == io.EOF {
		return
	} else {
		common.ExitIfError(err)
	}
}

// Close closes the File, rendering it unusable for I/O. On files that support SetDeadline,
// any pending I/O operations will be canceled and return immediately with an error.
// Close will return an error if it has already been called.
func (reader *SimpleReader) Close() {
	if reader != nil {
		err := reader.File.Close()
		common.ExitIfError(err)
	}
}
