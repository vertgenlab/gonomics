package sam

import (
	"os"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/test.sam"},
}

// TODO: better read test
func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		Read(test.filename)
	}
}

// TODO: better rw test
func TestReadAndWrite(t *testing.T) {
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		actual, header := Read(test.filename)
		Write(tempFile, actual, header)
		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}
