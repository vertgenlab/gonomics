package fileio

import (
	"bytes"
	"io"
	"os"
	"strings"
	"testing"

	"github.com/klauspost/pgzip"
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

func TestMustCreateEmptyFilenameError(t *testing.T) {
	if os.Getenv("TEST_FATAL") != "1" {
		return
	}
	// Will cause log.Fatal in this test context
	MustCreate("")
}

func TestStdinSimulation(t *testing.T) {
	testCases := []string{
		"Gonomics fileio stdin mock",
		"",                 // Empty string
		"\nHello\nWorld\n", // Multiple lines
		"Secret input with special characters !@#$%^&*()_+",
	}

	for _, input := range testCases {
		mockStdin, err := os.CreateTemp("", "mockStdin")
		if err != nil {
			t.Fatalf("Error: Failed to create temporary file: %v", err)
		}
		defer mockStdin.Close()
		if _, err := mockStdin.WriteString(input); err != nil {
			t.Fatalf("Error: Failed to write simulated input: %v", err)
		}
		if _, err := mockStdin.Seek(0, io.SeekStart); err != nil {
			t.Fatalf("Error: Failed to seek to the start: %v", err)
		}
		os.Stdin = mockStdin

		file := MustOpen("/dev/stdin")
		defer file.Close()

		var reader strings.Builder

		if _, err := io.Copy(&reader, file); err != nil {
			t.Fatalf("Error: Failed to read from simulated stdin: %v", err)
		}
		if reader.String() != input {
			t.Errorf("Error: Mismatch for input %q. Expected: %q, got %q", input, input, reader.String())
		}
		MustRemove(mockStdin.Name())
	}
}

func TestStdinGzip(t *testing.T) {
	testCases := []string{
		"fileio stdin with gzip",
		"stdin.gz",
	}

	for _, input := range testCases {
		// Create a temporary file to simulate stdin
		mockStdin, err := os.CreateTemp("", "mockStdin.gz")
		if err != nil {
			t.Fatalf("Error: Failed to create temporary file: %v", err)
		}
		defer os.Remove(mockStdin.Name())
		defer mockStdin.Close()

		// Gzip the input data
		var gzippedInput bytes.Buffer
		gzWriter := pgzip.NewWriter(&gzippedInput)
		if _, err := gzWriter.Write([]byte(input)); err != nil {
			t.Fatalf("Error: Failed to gzip input: %v", err)
		}
		if err := gzWriter.Close(); err != nil {
			t.Fatalf("Error: Failed to close gzip writer: %v", err)
		}

		// Write the gzipped input to the mock stdin
		if _, err := io.Copy(mockStdin, &gzippedInput); err != nil {
			t.Fatalf("Error: Failed to write simulated input: %v", err)
		}

		// Seek to the beginning of the mock stdin
		if _, err := mockStdin.Seek(0, io.SeekStart); err != nil {
			t.Fatalf("Error: Failed to seek to the start: %v", err)
		}

		// Set os.Stdin to the mock stdin
		os.Stdin = mockStdin

		// Call the function being tested
		file := EasyOpen("/dev/stdin")
		defer file.Close()

		// Read from the file and compare with the original input
		var reader strings.Builder
		if _, err := io.Copy(&reader, file); err != nil {
			t.Fatalf("Error: Failed to read from simulated stdin: %v", err)
		}
		if reader.String() != input {
			t.Errorf("Error: Mismatch for input %q. Expected: %q, got %q", input, input, reader.String())
		}
	}
}
