package fileio

import (
	"io"
	"os"
	"strings"
	"testing"
)

var testfile string = "testdata/smallTest"
var line1 string = "#shhhh this line is a secret"
var line2 string = "Hello World"
var line3 string = "I am a gopher"

func TestNextLine(t *testing.T) {
	file := EasyOpen(testfile)
	l1, done := NextLine(file.BuffReader)
	if l1 != line1 || done {
		t.Errorf("problem reading next line")
	}
	l2, done := NextLine(file.BuffReader)
	if l2 != line2 || done {
		t.Errorf("problem reading next line")
	}
	l3, done := NextLine(file.BuffReader)
	if l3 != line3 || done {
		t.Errorf("problem reading next line")
	}
	l4, done := NextLine(file.BuffReader)
	if l4 != "" || !done {
		t.Errorf("problem reading next line")
	}
}

func TestNextRealLine(t *testing.T) {
	file := EasyOpen(testfile)
	l1, done := NextRealLine(file.BuffReader)
	if l1 != line2 || done {
		t.Errorf("problem reading next line")
	}
	l2, done := NextRealLine(file.BuffReader)
	if l2 != line3 || done {
		t.Errorf("problem reading next line")
	}
	l3, done := NextRealLine(file.BuffReader)
	if l3 != "" || !done {
		t.Errorf("problem reading next line")
	}
}

func TestEqual(t *testing.T) {
	if !AreEqual(testfile, testfile) || !AreEqualIgnoreComments(testfile, testfile) {
		t.Errorf("problem with equal")
	}
}

func TestFileioStdin(t *testing.T) {
	input := "Gonomics!"

	// Create a temporary file to simulate input to stdin
	stdin, err := os.CreateTemp("", "stdinTest")
	if err != nil {
		t.Fatalf("Error: Failed to create temporary file: %v", err)
	}

	if _, err := stdin.WriteString(input); err != nil {
		t.Fatalf("Error: Failed to write simulated input to temp file: %v", err)
	}

	// Seek back to the start of the file so it can be read
	if _, err := stdin.Seek(0, io.SeekStart); err != nil {
		t.Fatalf("Error: Failed to seek to the start of the temp file: %v", err)
	}
	// Redirect os.Stdin to the temporary file
	os.Stdin = stdin

	// Call MustOpen with "stdin" and read from os.Stdin
	file := MustOpen("/dev/stdin")
	defer file.Close()

	// Read the data returned by MustOpen
	var readData strings.Builder
	if _, err := io.Copy(&readData, file); err != nil {
		t.Fatalf("Error: Failed to read from simulated stdin: %v", err)
	}

	// Verify the data matches what was written to the simulated stdin
	if readData.String() != input {
		t.Errorf("Error: Expected %q, got %q", input, readData.String())
	}
	os.Remove(stdin.Name())
}
