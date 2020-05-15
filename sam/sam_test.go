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

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		_, err := Read(test.filename)
		if err != nil {
			t.Errorf("Reading %s gave an error..", test.filename)
		}
	}
}

func TestReadAndWrite(t *testing.T) {
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		actual, err := Read(test.filename)
		if err != nil {
			t.Errorf("Error writing %s as a temp fasta file", tempFile)
		}
		Write(tempFile, actual)
		if err != nil {
			t.Errorf("Reading %s gave an error", test.filename)
		}
		/*if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
		*/
		err = os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}
