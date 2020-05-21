package axt

import (
	"os"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/chrM_gasacu1.axt"},
}

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		_ = Read(test.filename)
	}
}

func TestWriteAndRead(t *testing.T) {
	var actual []*Axt
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		actual = Read(test.filename)
		Write(tempFile, actual)
		if !AllEqual(Read(tempFile), Read("testdata/chrM_gasacu1.axt")) {
			t.Errorf("Axt files are not the same")
		}

		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}
