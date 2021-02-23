package maf

import (
	//"os"
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	//{"testdata/chr7.test.maf"},
	{"testdata/chr22.test.maf"},
}

func TestReadWrite(t *testing.T) {
	var actual []*Maf
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		actual = Read(test.filename)
		Write(tempFile, actual)
		if !fileio.AreEqualIgnoreComments(test.filename, tempFile) {
			t.Errorf("File has changed after reading and writing.")
		}
		fileio.EasyRemove(tempFile)
	}
}
