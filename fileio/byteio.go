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

// ByteReader provides enhanced byte reading from files, with transparent support for pgzip decompression
type ByteReader struct {
	*bufio.Reader
	File         *os.File
	internalGzip *pgzip.Reader // Pgzip reader for gzipped files.
	line         []byte
	Buffer       *bytes.Buffer // Internal line buffer.
}

// Read reads data into p and is a method required to implement the io.Reader interface.
func (reader *ByteReader) Read(b []byte) (n int, err error) {
	return reader.Read(b)
}

// NewByteReader creates a ByteReader for reading from the specified file, automatically handling gzip decompression
func NewByteReader(filename string) *ByteReader {
	var err error
	var answer = ByteReader{
		File:   MustOpen(filename),
		Buffer: &bytes.Buffer{},
	}
	switch true {
	case strings.HasSuffix(filename, ".gz"):
		answer.internalGzip, err = pgzip.NewReader(answer.File)
		exception.PanicOnErr(err)
		answer.Reader = bufio.NewReader(answer.internalGzip)
	default:
		answer.internalGzip = nil
		answer.Reader = bufio.NewReader(answer.File)
	}
	return &answer
}

// ReadLine reads a single line from the file into bufio.Reader, handling line breaks.
func ReadLine(reader *ByteReader) (*bytes.Buffer, bool) {
	var err error
	reader.line, err = reader.Reader.ReadBytes('\n') // Assuming Reader is a bufio.Reader or similar.
	if err != nil && err != io.EOF {
		// Handle non-EOF errors directly. For EOF, we'll handle it after checking the line content.
		log.Fatalf("Error reading line: %v", err)
	}
	// If line is not empty or err is EOF (might still have data to return)
	if len(reader.line) > 0 || err == io.EOF {
		// Check for and trim a trailing newline character, if present.
		if len(reader.line) > 0 && reader.line[len(reader.line)-1] == '\n' {
			reader.line = reader.line[:len(reader.line)-1] // Trim \n
			// Further check for \r\n (Windows) reader.line endings if necessary:
			if len(reader.line) > 0 && reader.line[len(reader.line)-1] == '\r' {
				reader.line = reader.line[:len(reader.line)-1] // Trim \r
			}
		}
		// Write the line to the buffer.
		reader.Buffer.Reset()
		_, err = reader.Buffer.Write(reader.line)
		if err != nil {
			exception.PanicOnErr(err)
		}
		// Return the buffer and a flag indicating whether we should continue reading.
		// If err is io.EOF, this might be the last line.
		return reader.Buffer, err != io.EOF
	}
	// If we got here because of an EOF with no data, return nil and false.
	return nil, false
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

// Close closes the ByteReader, releasing any resources associated with the file and the gzip reader, if present.
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
