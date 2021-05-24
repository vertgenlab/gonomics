package wig

import (
	"os"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/in_test.wig"},
}

func TestWriteAndRead(t *testing.T) {
	var actual []Wig
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		actual = Read(test.filename)
		Write(tempFile, actual)

		if !AllEqual(Read(tempFile), Read("testdata/in_test.wig")) {
			t.Errorf("Read and write Wig files were not the same")
		}
		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}
