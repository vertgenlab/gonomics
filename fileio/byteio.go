package fileio

import (
	"bufio"
	"bytes"
	"io"
	"log"
	"os/exec"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/exception"
)

const (
	defaultBufSize = 4096
)

// ByteReader implements the io.Reader interface by providing the Read(b []byte) method.
// It wraps a *bufio.Reader and optionally a GunzipReader for reading gzipped files.
// The primary advantage of ByteReader is that it allows data to be read into a shared bytes.Buffer
// instead of allocating individual strings.
type ByteReader struct {
	*bufio.Reader
	internalGzip *GunzipReader
	line         []byte
	Buffer       *bytes.Buffer
	close        func() error
}

// GunzipReader implements the io.ReadCloser interface and is used to read gzipped files.
// It wraps an io.ReadCloser and a command for executing the gunzip command.
type GunzipReader struct {
	io.ReadCloser
	Cmd *exec.Cmd
}

// Read reads data into p and implements the io.Reader interface for ByteReader.
func (reader *ByteReader) Read(b []byte) (n int, err error) {
	return reader.Reader.Read(b)
}

// Read reads data into p and implements the io.Reader interface for GunzipReader.
func (reader *GunzipReader) Read(b []byte) (n int, err error) {
	return reader.ReadCloser.Read(b)
}

// NewByteReader creates a new ByteReader for the given filename.
// It automatically detects if the file is gzipped and handles it accordingly.
func NewByteReader(filename string) *ByteReader {
	var answer ByteReader
	answer.Buffer = &bytes.Buffer{}
	switch {
	case strings.HasSuffix(filename, ".gz"):
		answer.internalGzip = NewGunzipReader(filename)
		answer.Reader = bufio.NewReader(answer.internalGzip)
		answer.close = answer.internalGzip.Close
	default:
		file := MustOpen(filename)
		answer.Reader = bufio.NewReader(file)
		answer.close = file.Close
	}
	return &answer
}

// NewGunzipReader creates a new GunzipReader for the given filename.
// It executes the gunzip command to read the gzipped file.
func NewGunzipReader(filename string) *GunzipReader {
	var answer GunzipReader
	answer.Cmd = exec.Command("gunzip", "-c", filename)
	stdout, err := answer.Cmd.StdoutPipe()
	exception.PanicOnErr(err)
	answer.ReadCloser = stdout
	exception.PanicOnErr(answer.Cmd.Start())
	return &answer
}

// ReadLine reads one line from the ByteReader and returns it as a bytes.Buffer.
// It also returns a bool indicating whether there are more lines to read.
// The returned buffer points to the internal slice belonging to the reader, so it should be copied if needed.
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
			reader.line = recursiveRead(reader)
			return bytesToBuffer(reader), false
		} else {
			catchErrThrowEOF(err)
		}
	}
	return nil, true
}

// recursiveRead is a helper function to deal with very long lines and avoid allocating too much memory upfront.
// It recursively reads the next bytes until reaching the end of line character.
func recursiveRead(reader *ByteReader) []byte {
	_, err := reader.Buffer.Write(reader.line)
	exception.PanicOnErr(err)
	reader.line, err = reader.ReadSlice('\n')
	if err == nil {
		return reader.line
	}
	if err == bufio.ErrBufferFull {
		_, err = reader.Buffer.Write(reader.line)
		exception.PanicOnErr(err)
		return recursiveRead(reader)
	}
	exception.PanicOnErr(err)
	return reader.line
}

// catchErrThrowEOF silently handles and throws the EOF error, and logs and exits on any other error.
func catchErrThrowEOF(err error) {
	if err == io.EOF {
		return
	} else {
		exception.PanicOnErr(err)
	}
}

// bytesToBuffer converts the internal byte slice of the ByteReader to a bytes.Buffer and returns it.
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

// Close closes the ByteReader, rendering it unusable for I/O.
func (br *ByteReader) Close() error {
	if br.internalGzip != nil {
		exception.PanicOnErr(br.internalGzip.Close())
	} else {
		exception.PanicOnErr(br.close())
	}
	return nil
}

// Close closes the GunzipReader, rendering it unusable for I/O.
func (unzip *GunzipReader) Close() error {
	exception.PanicOnErr(unzip.ReadCloser.Close())
	exception.PanicOnErr(unzip.Cmd.Wait())
	return nil
}

// StringToIntSlice converts a comma-separated string to a slice of integers.
// It handles formats like PSL and genePred that may have a trailing comma.
func StringToIntSlice(line string) []int {
	work := strings.Split(line, ",")
	sliceSize := len(work)
	if line[len(line)-1] == ',' {
		sliceSize--
	}
	answer := make([]int, sliceSize)
	for i := 0; i < sliceSize; i++ {
		num, err := strconv.Atoi(work[i])
		exception.PanicOnErr(err)
		answer[i] = num
	}
	return answer
}

// IntSliceToString converts a slice of integers to a comma-separated string.
// The string will include a trailing comma to satisfy UCSC's anomalies.
func IntSliceToString(nums []int) string {
	ans := strings.Builder{}
	ans.Grow(2 * len(nums))
	for i := 0; i < len(nums); i++ {
		ans.WriteString(strconv.Itoa(nums[i]))
		ans.WriteByte(',')
	}
	return ans.String()
}

// IntToString converts an integer to a string.
func IntToString(i int) string {
	return strconv.Itoa(i)
}
