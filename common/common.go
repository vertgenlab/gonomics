package common

import (
	"bufio"
	"log"
	"strings"
	"os"
)

func ExitIfError(err error) {
	if err != nil {
		log.Fatal(err)
	}
}

func Exit(message string) {
	log.Fatal(message)
}

func MustOpen(filename string) *os.File {
	file, err := os.Open(filename)
	ExitIfError(err)
	return file
}

// returns the next line of the file (might be a comment line)
func NextLine(reader *bufio.Reader) (string, error) {
	var line string
	var err error
	for line, err = reader.ReadString('\n'); err == nil; line, err = reader.ReadString('\n') {
	}
	line = strings.TrimSuffix(line, "\n")
	return line, err
}

// returns the next line of the file that is not a comment line
func NextRealLine(reader *bufio.Reader) (string, error) {
	var line string
	var err error
	for line, err = reader.ReadString('\n'); err == nil && strings.HasPrefix(line, "#"); line, err = reader.ReadString('\n') {
	}
	line = strings.TrimSuffix(line, "\n")
	return line, err
}

func Max(a int, b int) int {
	if a >= b {
		return a
	} else {
		return b
	}
}

func Min(a int, b int) int {
	if a <= b {
		return a
	} else {
		return b
	}
}

func TripleMax(a int, b int, c int) int {
	if a >= b && a >= c {
		return a
	} else if b >= c {
		return b
	} else {
		return c
	}
}

func TripleMin(a int, b int, c int) int {
	if a <= b && a <= c {
		return a
	} else if b <= c {
		return b
	} else {
		return c
	}
}
