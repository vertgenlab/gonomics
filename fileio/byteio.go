package fileio

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"errors"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

const (
	defaultBufSize = 4096
)

// ByteReader implements the io.Reader interface by providing
// the Read(b []byte) method. The struct contains an embedded *bufio.Reader
// and a pointer to os.File for closure when reading is complete.
// The primary advantage of ByteReader over EasyReader is that data
// is read into a shared bytes.Buffer instead of allocating memory for a
// string as is done with EasyReader.
// ByteReader can be be used to efficiently parse files on the byte level
// which may be significantly faster and require less memory than using
// strings.Split on a string derived from EasyReader.
// The drawback is that ByteReader is not as easy to use as EasyReader
// and should generally be reserved for performance intensive tasks.
type ByteReader struct {
	*bufio.Reader
	File         *os.File
	internalGzip *gzip.Reader
	line         []byte
	Buffer       *bytes.Buffer
}

// Read reads data into p and is a method required to implement the io.Reader interface.
// It returns the number of bytes read into p.
func (reader *ByteReader) Read(b []byte) (n int, err error) {
	return reader.Read(b)
}

// NewByteReader will process a given file and performs error handling if an error occurs.
// ByteReader will process gzipped files accordingly by performing a check on the suffix
// of the provided file.
func NewByteReader(filename string) *ByteReader {
	var err error
	file := MustOpen(filename)
	var answer ByteReader = ByteReader{
		File:   file,
		Buffer: &bytes.Buffer{},
	}
	switch true {
	case strings.HasSuffix(filename, ".gz"):
		answer.internalGzip, err = gzip.NewReader(file)
		exception.PanicOnErr(err)
		answer.Reader = bufio.NewReader(answer.internalGzip)
	default:
		answer.Reader = bufio.NewReader(file)
	}
	return &answer
}

// ReadLine will return a bytes.Buffer pointing to the internal slice of bytes. Provided this function is called within a loop,
// the function will read one line at a time, and return bool to continue reading. Important to note the buffer return points to
// the internal slice belonging to the reader, meaning the slice will be overridden if the data is not copied.
func ReadLine(reader *ByteReader) (*bytes.Buffer, bool) {
	var err error
	reader.line, err = reader.ReadSlice('\n')
	reader.Buffer.Reset()
	if err == nil {
		if reader.line[len(reader.line)-1] == '\n' {
			return bytesToBuffer(reader), false
		} else {
			log.Panicf("Error: end of line did not end with an end of line character...\n")
		}
	} else {
		if err == bufio.ErrBufferFull {
			reader.line = readMore(reader)
			return bytesToBuffer(reader), false
		} else {
			CatchErrThrowEOF(err)
		}
	}
	return nil, true
}

// readMore is a private helper function to deal with very long lines to
// avoid alocating too much memory upfront and only resize the size of the buffer
// only when necessary.
func readMore(reader *ByteReader) []byte {
	_, err := reader.Buffer.Write(reader.line)
	exception.PanicOnErr(err)
	reader.line, err = reader.ReadSlice('\n')
	if err == nil {
		return reader.line
	}
	if err == bufio.ErrBufferFull {
		_, err = reader.Buffer.Write(reader.line)
		exception.PanicOnErr(err)
		// recursive call to read next bytes until reaching end of line character
		return readMore(reader)
	}
	exception.PanicOnErr(err)
	return reader.line
}

// CatchErrThrowEOF will silently handles and throws the EOF error and will log and exit any other errors.
func CatchErrThrowEOF(err error) {
	if err == io.EOF {
		return
	} else {
		exception.PanicOnErr(err)
	}
}

// bytesToBuffer will parse []byte and return a pointer to the same underlying bytes.Buffer
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

// Close closes the File, rendering it unusable for I/O. On files that support SetDeadline,
// any pending I/O operations will be canceled and return immediately with an error.
// Close will return an error if it has already been called.
func (br *ByteReader) Close() error {
	var gzErr, fileErr error
	if br.internalGzip != nil {
		gzErr = br.internalGzip.Close()
	}
	if br.File != nil {
		fileErr = br.File.Close()
	} else {
		return errors.New("no file found")
	}

	switch { // Handle error returns. Priority is gzErr > fileErr
	case gzErr != nil:
		return gzErr

	case fileErr != nil:
		log.Println("WARNING: attempted to close file, but file already closed")
		return nil

	default:
		return nil
	}
}

// StringToIntSlice will process a row of data separated by commas, convert the slice into a slice of type int.
// PSL and genePred formats have a trailing comma we need to account for and the check at the beginning will adjust
// the length of the working slice.
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

// IntListToString will process a slice of type int as an input and return a each value separated by a comma as a string.
// Important Note: string will include a trailing comma to satisfy UCSC's anomalies.
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
