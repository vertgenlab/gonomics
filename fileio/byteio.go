package fileio

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/klauspost/pgzip"
	"github.com/vertgenlab/gonomics/exception"
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
	internalGzip *pgzip.Reader // Use pgzip.Reader here
	line         []byte
	Buffer       *bytes.Buffer
}

// Read reads data into p and is a method required to implement the io.Reader interface.
// It returns the number of bytes read into p.
func (reader *ByteReader) Read(b []byte) (n int, err error) {
	if reader.internalGzip != nil {
		return reader.internalGzip.Read(b)
	}
	return reader.Reader.Read(b)
}

// NewByteReader will process a given file and performs error handling if an error occurs.
// ByteReader will process gzipped files accordingly by performing a check on the suffix
// of the provided file.
func NewByteReader(filename string) *ByteReader {
	file := MustOpen(filename)

	reader := &ByteReader{
		File:   file,
		Buffer: &bytes.Buffer{},
	}

	if strings.HasSuffix(filename, ".gz") {
		var err error
		reader.internalGzip, err = pgzip.NewReader(file)
		exception.PanicOnErr(err)

		reader.Reader = bufio.NewReader(reader.internalGzip)
	} else {
		reader.Reader = bufio.NewReader(file)
	}
	return reader
}

// ReadLine will return a bytes.Buffer pointing to the internal slice of bytes. Provided this function is called within a loop,
// the function will read one line at a time, and return bool to continue reading. Important to note the buffer return points to
// the internal slice belonging to the reader, meaning the slice will be overridden if the data is not copied.
func ReadLine(reader *ByteReader) (*bytes.Buffer, bool) {
	reader.Buffer.Reset() // Reset buffer for new line reading
	var err error
	for reader.line, err = reader.ReadSlice('\n'); err != nil || err != bufio.ErrBufferFull; reader.line, err = reader.ReadSlice('\n') {
		if err == io.EOF {
			return reader.Buffer, true // End of file reached
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
			break // Complete line read, exit loop
		} else {
			_, err = reader.Buffer.Write(reader.line)
			exception.PanicOnErr(err)
			if err == bufio.ErrBufferFull {
				continue // Buffer was full, continue reading the line
			}
		}
	}
	return reader.Buffer, false
}

// CatchErrThrowEOF will silently handles and throws the EOF error and will log and exit any other errors.
func CatchErrThrowEOF(err error) {
	if err == io.EOF {
		return
	} else {
		exception.PanicOnErr(err)
	}
}

// bytesToBuffer will parse []byte and return a pointer to the same underlying bytes.Buffer.
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
