package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var BedSubsetMatrixTests = []struct {
	UnionFile                    string
	FileListFile                 string
	Fraction                     bool
	MarkMultipleOverlaps         string
	OutFile                      string
	ExpectedFile                 string
	ExpectedMultipleOverlapsFile string
}{
	{"testdata/union.bed", "testdata/files.list", false, "", "testdata/out1.txt", "testdata/expected.txt", ""},
	{"testdata/union.bed", "testdata/files.fraction1.list", true, "", "testdata/out2.txt", "testdata/expected.fraction1.txt", ""},
	{"testdata/union.bed", "testdata/files.fraction2.list", true, "", "testdata/out3.txt", "testdata/expected.fraction2.txt", ""},
	{"testdata/union.bed", "testdata/files.fraction3.list", false, "", "testdata/out4.txt", "testdata/expected.noFraction2.txt", ""},
	{"testdata/union.bed", "testdata/files.fraction3.list", true, "", "testdata/out5.txt", "testdata/expected.fraction3.txt", ""},
	{"testdata/union.bed", "testdata/files.fraction3.list", true, "testdata/out5.multipleOverlaps.txt", "testdata/out5.txt", "testdata/expected.fraction3.txt", "testdata/expected.fraction3.multipleOverlaps.txt"},
}

func TestBedSubsetMatrix(t *testing.T) {
	var err error

	for _, v := range BedSubsetMatrixTests {
		intervalSubsetMatrix(v.UnionFile, v.FileListFile, v.OutFile, v.Fraction, v.MarkMultipleOverlaps)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in bedSubsetMatrix. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
		if v.MarkMultipleOverlaps != "" {
			if !fileio.AreEqual(v.MarkMultipleOverlaps, v.ExpectedMultipleOverlapsFile) {
				t.Errorf("Error in bedSubsetMatrix -markMultipleOverlaps. Multiple overlaps output file was not as expected.")
			} else {
				err = os.Remove(v.MarkMultipleOverlaps)
				exception.PanicOnErr(err)
			}
		}
	}
}
