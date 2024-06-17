package fileio

import (
	"bufio"
	"bytes"
	"encoding/binary"
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

// ByteReader wraps bufio.Reader for efficient byte-level file parsing.
type ByteReader struct {
	*bufio.Reader
	File         *os.File
	internalGzip *pgzip.Reader
	line         []byte
	Buffer       *bytes.Buffer
}

// Read implements io.Reader, reading data into b from the file or gzip stream.
func (reader *ByteReader) Read(b []byte) (n int, err error) {
	if reader.internalGzip != nil {
		return reader.internalGzip.Read(b)
	}
	return reader.Reader.Read(b)
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

// Close closes the internal File and gzip Reader if present.
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

// DecodeLittleEndianBinaryField reads little endian binary data from an input *EasyReader to a variable.
func DecodeLittleEndianBinaryField(file io.Reader, data any) {
	err := binary.Read(file, binary.LittleEndian, data)
	exception.PanicOnErr(err)
}
