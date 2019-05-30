package common

import (
	"bufio"
	"io"
	"os"
	"strings"
)

func MustOpen(filename string) *os.File {
	file, err := os.Open(filename)
	ExitIfError(err)
	return file
}

// returns the next line of the file (might be a comment line)
// returns true if the file is done
func NextLine(reader *bufio.Reader) (string, bool) {
	var line string
	var err error
	line, err = reader.ReadString('\n')
	if err != nil && err != io.EOF {
		ExitIfError(err)
	}
	line = strings.TrimSuffix(line, "\n")
	if err == io.EOF {
		return line, true
	}
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
		ExitIfError(err)
	}
	line = strings.TrimSuffix(line, "\n")
	if err == io.EOF {
		return line, true
	}
	return line, false
}
