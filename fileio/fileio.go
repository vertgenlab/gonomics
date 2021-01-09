package fileio

import (
	"bufio"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"io"
	"log"
	"os"
	"strings"
)

func MustCreate(filename string) *os.File {
	file, err := os.Create(filename)
	common.ExitIfError(err)
	return file
}

func MustOpen(filename string) *os.File {
	file, err := os.Open(filename)
	common.ExitIfError(err)
	return file
}

// returns the next line of the file (might be a comment line)
// returns true if the file is done
func NextLine(reader *bufio.Reader) (string, bool) {
	var line string
	var err error
	line, err = reader.ReadString('\n')
	if err != nil && err != io.EOF {
		common.ExitIfError(err)
	}
	if err == io.EOF {
		if line != "" {
			log.Fatalf("Error: last line of file didn't end with a newline character: %s\n", line)
		} else {
			return "", true
		}
	}
	line = strings.TrimSuffix(line, "\n")
	return line, false
}

// returns the next line of the file that is not a comment line
// returns true if the file is done
func NextRealLine(reader *bufio.Reader) (string, bool) {
	var line string
	var err error
	for line, err = reader.ReadString('\n'); err == nil && strings.HasPrefix(line, "#"); line, err = reader.ReadString('\n') {
	}
	if err != nil && err != io.EOF {
		common.ExitIfError(err)
	}
	if err == io.EOF {
		if line != "" {
			log.Fatalf("Error: last line of file didn't end with a newline character: %s\n", line)
		} else {
			return "", true
		}
	}
	line = strings.TrimSuffix(line, "\n")
	return line, false
}

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

func AreEqualIgnoreComments(a string, b string) bool {
	return equal(a, b, false)
}

func AreEqual(a string, b string) bool {
	return equal(a, b, true)
}

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
