package fileio

import (
	"bufio"
	"bytes"
	"fmt"
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
// - The primary advantage of ByteReader is data read into a shared bytes.Buffer rather than allocating memory.
// - The drawback is that ByteReader is not as easy to use as EasyReader nd should be reserved for performance intensive tasks.
// The struct contains:
// - *bufio.Reader, embedded
// - GunzipReader which implements io.ReadCloser, a builtin io.Reader and io.Closer interface.
// - []byte
// - *bytes.Buffer
// - close function from os.File
type ByteReader struct {
	*bufio.Reader
	internalGzip *GunzipReader
	line         []byte
	Buffer       *bytes.Buffer
	close        func() error
}

// GunzipReader is implements the io.ReadCloser interface, which is a golang buildin interface that groups the basic Read and Close methods.
type GunzipReader struct {
	io.ReadCloser
	Cmd *exec.Cmd
}

// Read reads data into p and is a method required by ByteReader implement the io.Reader interface.
func (reader *ByteReader) Read(b []byte) (n int, err error) {
	return reader.Read(b)
}

// Read reads data into p and is a method required by GunzipReader implement the io.Reader interface.
func (reader GunzipReader) Read(b []byte) (n int, err error) {
	return reader.ReadCloser.Read(b)
}

// NewByteReader will process a given file and performs error handling if an error occurs and process gzipped files accordingly by suffix of the provided file.
func NewByteReader(filename string) *ByteReader {
	var answer ByteReader = ByteReader{
		Buffer: &bytes.Buffer{},
		line:   make([]byte, defaultBufSize),
	}
	switch true {
	case strings.HasSuffix(filename, ".gz"):
		answer.internalGzip = NewGunzipReader(filename)
		answer.Reader = bufio.NewReader(answer.internalGzip)
		answer.close = answer.internalGzip.ReadCloser.Close
	default:
		file := MustOpen(filename)
		answer.Reader = bufio.NewReader(file)
		answer.close = file.Close
	}
	return &answer
}

// NewGunzipReader takes a file name and returns an implementation of io.Reader.
func NewGunzipReader(filename string) *GunzipReader {
	var answer GunzipReader = GunzipReader{
		Cmd: exec.Command("gunzip", "-c", filename),
	}
	stdout, err := answer.Cmd.StdoutPipe()
	exception.PanicOnErr(err)
	answer.ReadCloser = stdout
	exception.PanicOnErr(answer.Cmd.Start())
	return &answer
}

// ReadLine will return a bytes.Buffer pointing to the internal slice of bytes.
// Important to note:
// - Provided this function is called within a loop, the function will read one line at a time, and return bool to continue reading.
// - the buffer return points to the internal slice belonging to the reader, meaning the slice will be overridden if the data is not copied.
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
			CatchErrThrowEOF(err)
		}
	}
	return nil, true
}

// recursiveRead is a private helper function to deal with very long lines to in an effort to avoid alocating too much memory upfront and only resize the size of the buffer when necessary.
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
		// recursive call to read next bytes until reaching end of line character
		return recursiveRead(reader)
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

// Close is called by ByteReader to closes files and render unusable for I/O. On files that support SetDeadline, any pending I/O operations will be canceled and return immediately with an error.
func (br *ByteReader) Close() error {
	if br.internalGzip != nil {
		exception.PanicOnErr(br.internalGzip.Close())
	} else {
		exception.PanicOnErr(br.close())
	}
	return nil
}

// Close is called by GunzipReader to closes files and render unusable for I/O. On files that support SetDeadline, any pending I/O operations will be canceled and return immediately with an error.
func (unzip *GunzipReader) Close() error {
	exception.PanicOnErr(unzip.ReadCloser.Close())
	exception.PanicOnErr(unzip.Cmd.Wait())
	return nil
}

// StringToIntSlice will process a row of data separated by commas, convert the slice into a slice of type int. PSL and genePred formats have a trailing comma we need to account for and the check at the beginning will adjust the length of the working slice.
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
// - Important Note: string will include a trailing comma to satisfy UCSC's anomalies.
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
