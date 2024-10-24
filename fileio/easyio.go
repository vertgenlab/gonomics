package fileio

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/klauspost/pgzip"
	"github.com/vertgenlab/gonomics/exception"
)

// EasyReader provides a simplified wrapper of the builtin golang io functions.
// Will silently handle reading gziped files. EasyReader is a valid io.Reader.
type EasyReader struct {
	File         *os.File
	internalGzip *pgzip.Reader
	BuffReader   *bufio.Reader
}

// EasyWriter provides a simplified wrapper of the builtin golang io functions.
// Will silently gzip output files when the input filename ends with '.gz'.
// EasyWriter is a valid io.Writer.
type EasyWriter struct {
	file         *os.File
	internalBuff *bufio.Writer
	internalGzip *pgzip.Writer
}

// EasyOpen opens the input file. Panics if errors are encountered.
func EasyOpen(filename string) *EasyReader {
	if strings.Contains(filename, "http") {
		return EasyHttp(filename)
	}
	answer := EasyReader{
		File: MustOpen(filename),
	}
	switch {
	case IsGzip(answer.File):
		var err error
		answer.internalGzip, err = pgzip.NewReader(answer.File)
		exception.PanicOnErr(err)
		answer.BuffReader = bufio.NewReader(answer.internalGzip)
	default:
		answer.BuffReader = bufio.NewReader(answer.File)
	}
	return &answer
}

// EasyCreate creates a file with the input name. Panics if errors are encountered.
func EasyCreate(filename string) *EasyWriter {
	answer := EasyWriter{}

	switch {
	case strings.HasPrefix(filename, "stdout"):
		answer.file = os.Stdout
	case strings.HasPrefix(filename, "stderr"):
		answer.file = os.Stderr
	default:
		answer.file = MustCreate(filename)
	}

	answer.internalBuff = bufio.NewWriter(answer.file)

	if strings.HasSuffix(filename, ".gz") {
		answer.internalGzip = pgzip.NewWriter(answer.internalBuff)
	} else {
		answer.internalGzip = nil
	}
	return &answer
}

// EasyNextLine returns the next line of the input EasyReader.
// Returns true at EOF.
func EasyNextLine(file *EasyReader) (string, bool) {
	return NextLine(file.BuffReader)
}

// EasyNextRealLine returns the next line of the input EasyReader that does not begin with '#'.
// Returns true at EOF.
func EasyNextRealLine(file *EasyReader) (string, bool) {
	return NextRealLine(file.BuffReader)
}

// EasyPeekReal will advance a reader past any lines beginning with '#' and read the first n bytes without advancing the reader.
func EasyPeekReal(file *EasyReader, n int) ([]byte, error) {
	return PeekReal(file.BuffReader, n)
}

// EasyReadHeader will read any leading comments lines from the file and return them as a slice of strings
// with each element in the slice being a comment line.
func EasyReadHeader(file *EasyReader) ([]string, error) {
	return ReadHeader(file.BuffReader)
}

// EasyRemove deletes the input file.
func EasyRemove(filename string) {
	MustRemove(filename)
}

// Read inputs a file and returns each line in the file as a string.
func Read(filename string) []string {
	var answer []string
	var err error
	file := EasyOpen(filename)
	reader := bufio.NewReader(file)
	for line, doneReading := NextRealLine(reader); !doneReading; line, doneReading = NextRealLine(reader) {
		answer = append(answer, line)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	return answer
}

// WriteToFileHandle will write a string to an io.Writer and panic on
// any error.
func WriteToFileHandle(file io.Writer, rec string) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", rec)
	exception.PanicOnErr(err)
}

// Write writes a slice of strings to a file with a newline
// placed after each element of the slice.
func Write(filename string, records []string) {
	var err error
	file := EasyCreate(filename)

	for i := range records {
		WriteToFileHandle(file, records[i])
	}
	err = file.Close()
	exception.PanicOnErr(err)
}
