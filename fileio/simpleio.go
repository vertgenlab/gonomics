package fileio

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"io"
	"log"
	"os"
	"strings"
)

const (
	defaultBufSize = 4096
)

// SimpleReader implements the io.Reader interface by providing
// the Read(b []byte) method. The struct contains an embedded *bufio.Reader
// and a pointer to os.File for closeure when reading is complete.
type SimpleReader struct {
	*bufio.Reader
	file   *os.File
	line   []byte
	Buffer *bytes.Buffer
}

// Read reads data into p and is a method required to implement the io.Reader interface.
// It returns the number of bytes read into p.
func (reader *SimpleReader) Read(b []byte) (n int, err error) {
	return reader.Read(b)
}

// NewSimpleReader will process a given file and performs error handling if an error occurs.
// SimpleReader will prcoess gzipped files accordinging by performing a check on the suffix
// of the provided file.
func NewSimpleReader(filename string) *SimpleReader {
	var answer SimpleReader = SimpleReader{
		file:   MustOpen(filename),
		line:   make([]byte, defaultBufSize),
		Buffer: &bytes.Buffer{},
	}
	switch true {
	case strings.HasSuffix(filename, ".gz"):
		gzipReader, err := gzip.NewReader(answer.file)
		common.ExitIfError(err)
		answer.Reader = bufio.NewReader(gzipReader)
	default:
		answer.Reader = bufio.NewReader(answer.file)
	}
	return &answer
}

// ReadLine will return a bytes.Buffer pointing to the internal slice of bytes. Provided this function is called within a loop,
// the function will read one line at a time, and return bool to continue reading. Important to note the buffer return points to
// the internal slice belonging to the reader, meaning the slice will be overridden if the data is not copied. Please be aware the
// reader will call close on the file once the reader encounters EOF.
func ReadLine(reader *SimpleReader) (*bytes.Buffer, bool) {
	var err error
	reader.line = reader.line[:0]
	reader.line, err = reader.ReadSlice('\n')
	if err == nil {
		if reader.line[len(reader.line)-1] == '\n' {
			reader.Buffer.Reset()
			_, err = reader.Buffer.Write(reader.line[:len(reader.line)-1])
			common.ExitIfError(err)
			return reader.Buffer, false
		} else {
			log.Fatalf("Error: end of line did not end with an end of line character...\n")
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
		err := reader.file.Close()
		common.ExitIfError(err)
	}
}

// StringToIntSlice will process a column data separated by commas, convert the slice into a slice of type int.
func StringToIntSlice(column string) []int {
	work := strings.Split(column, ",")
	var answer []int = make([]int, len(work))
	for i := 0; i < len(work)-1; i++ {
		answer[i] = common.StringToInt(work[i])
	}
	return answer
}

// IntListToString will process a slice of type int as an input and return a each value separated by a comma as a string.
func IntSliceToString(nums []int) string {
	ans := strings.Builder{}
	ans.Grow(2 * len(nums))
	for i := 0; i < len(nums); i++ {
		ans.WriteString(fmt.Sprintf("%d", nums[i]))
		ans.WriteByte(',')
	}
	return ans.String()
}

// IntToString a function that converts a number of type int and return a string.
func IntToString(i int) string {
	return fmt.Sprintf("%d", i)
}
