package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var MultFaVisualizerTests = []struct {
	Infile         string
	Outfile        string
	ExpectedFile   string
	Start          int
	End            int
	NoMask         bool
	Linelength     int
	EndOfAlignment bool
}{
	{"testdata/test.fa", "testdata/tmp.txt", "testdata/expected.txt", 1, 500, false, 50, false},
	{"testdata/test.fa", "testdata/tmp.noMask.txt", "testdata/expected.noMask.txt", 1, 500, true, 50, false},
	{"testdata/test.fa", "testdata/tmp.lineLength.txt", "testdata/expected.lineLength.txt", 1, 500, false, 100, false},
	{"testdata/test.fa", "testdata/tmp.short.txt", "testdata/expected.short.txt", 350, 400, false, 50, false},
	{"testdata/test.fa", "testdata/tmp.realShort.txt", "testdata/expected.realShort.txt", 4, 9, false, 50, false},
	{"testdata/test.fa", "testdata/tmp.4ToEnd.txt", "testdata/expected.4ToEnd.txt", 4, 9, false, 50, true},
}

func TestMultFaVisualizer(t *testing.T) {
	var err error
	for _, v := range MultFaVisualizerTests {
		multFaVisualizer(v.Infile, v.Outfile, v.Start, v.End, v.NoMask, v.Linelength, v.EndOfAlignment)
		if !fileio.AreEqual(v.Outfile, v.ExpectedFile) {
			t.Errorf("Error in MultFaVisualizer. Output does not match expected.")
		} else {
			err = os.Remove(v.Outfile)
			exception.PanicOnErr(err)
		}
	}
}
