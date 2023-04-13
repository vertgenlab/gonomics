package fileio

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"log"
	"os"
	"strings"

	"github.com/vertgenlab/gonomics/exception"
)

// EasyReader provides a simplified wrapper of the builtin golang io functions.
// Will silently handle reading gziped files. EasyReader is a valid io.Reader.
type EasyReader struct {
	File         *os.File
	internalGzip *gzip.Reader
	BuffReader   *bufio.Reader
}

// EasyWriter provides a simplified wrapper of the builtin golang io functions.
// Will silently gzip output files when the input filename ends with '.gz'.
// EasyWriter is a valid io.Writer.
type EasyWriter struct {
	file         *os.File
	internalBuff *bufio.Writer
	internalGzip *gzip.Writer
}

// EasyOpen opens the input file. Panics if errors are encountered.
func EasyOpen(filename string) *EasyReader {
	if strings.Contains(filename, "http") {
		return EasyHttp(filename)
	}

	answer := EasyReader{}
	var hasMagicGzip bool
	var readerInput io.Reader

	if strings.HasPrefix(filename, "stdin") {
		// when reading stdin we will assume the input is gzipped
		// if the file begins with the two magic gzip bytes 1f8d.
		// If it does, append .gz to the filename so it is parsed
		// as gzip in the following switch case.
		answer.File = os.Stdin
		readerInput, hasMagicGzip = newStdinMagicReader(magicGzip)
		if hasMagicGzip {
			filename += ".gz"
		}
	} else {
		answer.File = MustOpen(filename)
		hasMagicGzip = IsGzip(answer.File)
		readerInput = answer.File
	}

	var err error
	switch {
	case strings.HasSuffix(filename, ".gz") && hasMagicGzip:
		answer.internalGzip, err = gzip.NewReader(readerInput)
		exception.PanicOnErr(err)
		answer.BuffReader = bufio.NewReader(answer.internalGzip)

	case strings.HasSuffix(filename, ".gz"):
		log.Fatalf("ERROR: input file '%s' has the .gz suffix, but is not a gzip file", filename)

	case hasMagicGzip:
		log.Printf("WARNING: The input file '%s' looks like it may be gzipped, "+
			"but does not have the .gz suffix. Processing as a non-gzip file. Add the .gz "+
			"suffix to process as a gzip file.", filename)
		fallthrough

	default:
		answer.BuffReader = bufio.NewReader(readerInput)
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
		answer.internalGzip = gzip.NewWriter(answer.internalBuff)
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

func WriteToFileHandle(file io.Writer, rec string) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", rec)
	exception.PanicOnErr(err)
}

// Write writes a file.
func Write(filename string, records []string) {
	var err error
	file := EasyCreate(filename)

	for i := range records {
		WriteToFileHandle(file, records[i])
	}
	err = file.Close()
	exception.PanicOnErr(err)
}
