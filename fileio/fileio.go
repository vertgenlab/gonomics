// Package fileio provides wrappers of the builtin golang Reader/Writer utilities for ease of use and automatic gzip handling.
package fileio

import (
	"bufio"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"log"
	"os"
	"strings"
)

// MustCreate creates a file with the input name. Panics if not possible.
func MustCreate(filename string) *os.File {
	file, err := os.Create(filename)
	exception.PanicOnErr(err)
	return file
}

// MustOpen opens the input file. Panics if not possible.
func MustOpen(filename string) *os.File {
	file, err := os.Open(filename)
	exception.PanicOnErr(err)
	return file
}

// NextLine returns the next line of the file (might be a comment line).
// Returns true if the file is done.
func NextLine(reader *bufio.Reader) (string, bool) {
	var line string
	var err error
	line, err = reader.ReadString('\n')
	if err != nil && err != io.EOF {
		exception.PanicOnErr(err)
	}
	if err == io.EOF {
		if line != "" {
			log.Panicf("Error: last line of file didn't end with a newline character: %s\n", line)
		} else {
			return "", true
		}
	}
	line = strings.TrimSuffix(line, "\n")
	line = strings.TrimSuffix(line, "\r")
	return line, false
}

// NextRealLine returns the next line of the file that is not a comment line.
// Returns true if the file is done.
func NextRealLine(reader *bufio.Reader) (string, bool) {
	var line string
	var err error
	for line, err = reader.ReadString('\n'); err == nil && strings.HasPrefix(line, "#"); line, err = reader.ReadString('\n') {
	}
	if err != nil && err != io.EOF {
		log.Panic()
	}
	if err == io.EOF {
		if line != "" {
			log.Panicf("Error: last line of file didn't end with a newline character: %s\n", line)
		} else {
			return "", true
		}
	}
	line = strings.TrimSuffix(line, "\n")
	line = strings.TrimSuffix(line, "\r") //data generated from Windows OS contains \r\n as a two byte new line character.
	//Here we trim off trailing carriage returns. Lines without carriage returns are unaffected.
	return line, false
}

// PeekReal will advance a reader past any lines beginning with '#' and read the first n bytes without advancing the reader.
func PeekReal(reader *bufio.Reader, n int) ([]byte, error) {
	var peek []byte
	var err error
	for peek, err = reader.Peek(1); err == nil && peek[0] == '#'; peek, err = reader.Peek(1) {
		_, err = reader.ReadBytes('\n') // advance reader past comment line
		if err != nil {
			return nil, err
		}
	}

	if err != nil {
		return nil, err
	} else {
		return reader.Peek(n)
	}
}

// equal returns true if two input files are identical
func equal(a string, b string, commentsMatter bool) bool {
	var fileADone, fileBDone = false, false
	var lineA, lineB string

	fA := MustOpen(a)
	defer fA.Close()
	fB := MustOpen(b)
	defer fB.Close()
	readerA := bufio.NewReader(fA)
	readerB := bufio.NewReader(fB)

	for !fileADone && !fileBDone {
		if commentsMatter {
			lineA, fileADone = NextLine(readerA)
			lineB, fileBDone = NextLine(readerB)
		} else {
			lineA, fileADone = NextRealLine(readerA)
			lineB, fileBDone = NextRealLine(readerB)
		}
		if lineA != lineB {
			fmt.Printf("diff\n%s\n%s\n", lineA, lineB)
			return false
		}
	}
	if !fileADone || !fileBDone {
		return false
	}
	return true
}

// AreEqualIgnoreComments returns true if input files are equal.
// Ignores lines beginning with #.
func AreEqualIgnoreComments(a string, b string) bool {
	return equal(a, b, false)
}

// AreEqual returns true if input files are equal.
func AreEqual(a string, b string) bool {
	return equal(a, b, true)
}

// Read inputs a file and returns each line in the file as a string.
func Read(filename string) []string {
	var answer []string
	file := MustOpen(filename)
	defer file.Close()
	reader := bufio.NewReader(file)
	for line, doneReading := NextRealLine(reader); !doneReading; line, doneReading = NextRealLine(reader) {
		answer = append(answer, line)
	}
	return answer
}

//ReadFileToSingleLineString reads in any file type and returns contents without any \n
func ReadFileToSingleLineString(filename string) string {
	var catInput string
	var line string
	var doneReading bool = false
	file := EasyOpen(filename)
	defer file.Close()

	for line, doneReading = EasyNextRealLine(file); !doneReading; line, doneReading = EasyNextRealLine(file) {
		catInput = catInput + line
	}
	return catInput
}
