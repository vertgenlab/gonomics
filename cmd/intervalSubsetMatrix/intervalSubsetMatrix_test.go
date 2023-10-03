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
	Fraction     bool
	OutFile      string
	ExpectedFile string
}{
	{"testdata/union.bed", "testdata/files.list", false, "testdata/out1.txt", "testdata/expected.txt"},
	{"testdata/union.bed", "testdata/files.fraction1.list", true, "testdata/out2.txt", "testdata/expected.fraction1.txt"},
	{"testdata/union.bed", "testdata/files.fraction2.list", true, "testdata/out3.txt", "testdata/expected.fraction2.txt"},
}

func TestBedSubsetMatrix(t *testing.T) {
	var err error

	for _, v := range BedSubsetMatrixTests {
		intervalSubsetMatrix(v.UnionFile, v.FileListFile, v.OutFile, v.Fraction)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in bedSubsetMatrix. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
