package fileio

import (
	"bytes"
	"io"
	"log"
	"os"
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
)

var testfile string = "testdata/smallTest"
var gzipped string = "testdata/words.txt.gz"
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
func TestNotEqual(t *testing.T) {
	var buf bytes.Buffer
	log.SetOutput(&buf)
	defer log.SetOutput(os.Stderr)

	expected := "diff\nHello World\nmock text gzip\n"

	if equal(testfile, gzipped, true) || equal(testfile, gzipped, false) {
		t.Errorf("Error: results should not equal")
	}

	if !strings.HasSuffix(buf.String(), expected) {
		t.Errorf("Expected warning: '%s', got '%s'", expected, buf.String())
	}
}

func TestEqualGzip(t *testing.T) {
	var unzip string = "testdata/words.txt"

	copyFile(gzipped, unzip)
	defer EasyRemove(unzip)

	if !equal(gzipped, gzipped, true) || !equal(unzip, gzipped, true) {
		t.Errorf("Error: problem with unzip vs. gzip equal")
	}
}

func TestMustCreateEmptyFilenameError(t *testing.T) {
	var buf bytes.Buffer
	log.SetOutput(&buf)
	defer log.SetOutput(os.Stderr)
	defer exception.RecoverPanicErr()

	MustCreate("") // Will cause log.Panic in this test context
	logOutput := buf.String()

	expectedMessage := "must write to a non-empty filename"
	if !strings.Contains(expectedMessage, logOutput) {
		t.Errorf("Error: Expected log message:\n%s\nNot found in actual log:\n%s\n", expectedMessage, logOutput)
	} else {
		t.Log("log.Fatalf called with expected message - Test Passed!\n")
	}
}

func TestStdin(t *testing.T) {
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
		defer MustRemove(mockStdin.Name())
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
	}
}
