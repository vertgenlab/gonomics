package fileio

import (
	"bufio"
	"bytes"
	"io"
	"log"
	"os/exec"
	"strings"

	gzip "github.com/klauspost/pgzip"
	"github.com/vertgenlab/gonomics/exception"
)

// GunzipReader is bytesio that implements the system's native gzip.
type GunzipReader struct {
	*bufio.Reader
	Gunzip *zipio
	line   []byte
	Buffer *bytes.Buffer
	close  func() error
}

type ZipWriter struct {
	*bufio.Writer
	Gzip   *gzip.Writer
	Buffer *bytes.Buffer
	close  func() error
}

// zipio uncompress text using native system's gzip which out performs golang implementation.
type zipio struct {
	io.ReadCloser
	Cmd *exec.Cmd
}

// newZipReader takes a file name and returns an implementation of io.Reader.
func newZipReader(filename string) (*zipio, error) {
	cmd := exec.Command("gunzip", "-c", filename)
	stdout, err := cmd.StdoutPipe()
	exception.PanicOnErr(err)
	err = cmd.Start()
	return &zipio{ReadCloser: stdout, Cmd: cmd}, err
}

// NewGunzipReader will process a given file and performs error handling if an error occurs.
func NewGunzipReader(filename string) *GunzipReader {
	var answer GunzipReader = GunzipReader{
		line:   make([]byte, defaultBufSize),
		Buffer: &bytes.Buffer{},
	}
	switch true {
	case strings.HasSuffix(filename, ".gz"):
		var err error
		answer.Gunzip, err = newZipReader(filename)
		exception.PanicOnErr(err)
		answer.close = answer.Gunzip.ReadCloser.Close
		answer.Reader = bufio.NewReader(answer.Gunzip)
	default:
		answer.Gunzip = nil
		file := MustOpen(filename)
		answer.Reader = bufio.NewReader(file)
		answer.close = file.Close
	}
	return &answer
}

func GunzipLine(reader *GunzipReader) (*bytes.Buffer, bool) {
	var err error
	reader.line, err = reader.ReadSlice('\n')
	reader.Buffer.Reset()
	if err == nil {
		if reader.line[len(reader.line)-1] == '\n' {
			return GunzipToBuffer(reader), false
		} else {
			log.Panicf("Error: end of line did not end with an end of line character...\n")
		}
	} else {
		if err == bufio.ErrBufferFull {
			reader.line = readNext(reader)
			return GunzipToBuffer(reader), false
		} else {
			CatchErrThrowEOF(err)
		}
	}
	return nil, true
}

func readNext(reader *GunzipReader) []byte {
	_, err := reader.Buffer.Write(reader.line)
	exception.PanicOnErr(err)
	reader.line, err = reader.ReadSlice('\n')
	if err == nil {
		return reader.line
	}
	if err == bufio.ErrBufferFull {
		_, err = reader.Buffer.Write(reader.line)
		exception.PanicOnErr(err)
		return readNext(reader)
	}
	exception.PanicOnErr(err)
	return reader.line
}

func GunzipToBuffer(reader *GunzipReader) *bytes.Buffer {
	var err error
	_, err = reader.Buffer.Write(bytes.TrimSpace(reader.line))
	exception.PanicOnErr(err)
	return reader.Buffer
}

func (br GunzipReader) Close() error {
	var gzErr, fileErr error
	if br.Gunzip != nil {
		exception.PanicOnErr(br.close())
	}
	switch {
	case gzErr != nil:
		return gzErr
	case fileErr != nil:
		log.Println("Warning: attempted to close file, but file already closed")
		return nil
	default:
		return nil
	}
}

func (unzip zipio) Close() {
	unzip.Close()
	unzip.Cmd.Wait()
}

func NewWriter(filename string) *ZipWriter {
	ans := ZipWriter{}
	file := MustCreate(filename)
	ans.Writer = bufio.NewWriter(file)
	ans.Buffer = &bytes.Buffer{}
	if strings.HasSuffix(filename, ".gz") {
		ans.Gzip = gzip.NewWriter(ans.Writer)
		ans.Gzip.SetConcurrency(100000, 10)
	} else {
		ans.Gzip = nil
	}
	ans.close = file.Close
	return &ans
}

// SimplyRun uses a new goroutine to run the function.
func SimplyRun(f func()) {
	go f()
}
