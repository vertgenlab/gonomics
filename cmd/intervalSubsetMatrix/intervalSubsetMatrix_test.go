package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var BedSubsetMatrixTests = []struct {
	UnionFile    string
	FileListFile string
	OutFile      string
	ExpectedFile string
}{
	{"testdata/union.bed", "testdata/files.list", "testdata/out.txt", "testdata/expected.txt"},
}

func TestBedSubsetMatrix(t *testing.T) {
	var err error

	for _, v := range BedSubsetMatrixTests {
		intervalSubsetMatrix(v.UnionFile, v.FileListFile, v.OutFile)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in bedSubsetMatrix. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
