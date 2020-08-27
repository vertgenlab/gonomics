package fileio

import (
	"bufio"
	"compress/gzip"
	"github.com/vertgenlab/gonomics/common"
	"io"
	"log"
	"os"
	"strings"
)

// SimpleReader implements the io.Reader interface by providing
// the Read(b []byte) method. The struct contains an embedded *bufio.Reader
// and a pointer to os.File for closeure when reading is complete.
type SimpleReader struct {
	*bufio.Reader
	File *os.File
}

// Read reads data into p. It returns the number of bytes read into p.
func (reader *SimpleReader) Read(b []byte) (n int, err error) {
	return reader.Read(b)
}

// NewSimpleReader will process a given file and performs error handling if an error occurs.
// SimpleReader will prcoess gzipped files accordinging by performing a check on the suffix
// of the provided file.
func NewSimpleReader(filename string) *SimpleReader {
	var answer SimpleReader = SimpleReader{
		File: MustOpen(filename),
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

// ReadLine will return a slice of bytes for each line and provided this function is called within a loop,
// will return bool to continue reading. Please be aware that the function will call close on the file once
// the reader encounters EOF.
func ReadLine(reader *SimpleReader) ([]byte, bool) {
	curr, err := reader.ReadBytes('\n')
	if err == nil {
		if curr[len(curr)-1] == '\n' {
			return curr[:len(curr)-1], false
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
