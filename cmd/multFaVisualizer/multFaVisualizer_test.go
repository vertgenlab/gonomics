package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var MultFaVisualizerTests = []struct {
	Infile       string
	Outfile      string
	ExpectedFile string
	Start        int
	End          int
	NoMask       bool
	Linelength   int
}{
	{"testdata/test.fa", "testdata/tmp.txt", "testdata/expected.txt", 1, 500, false, 50},
	{"testdata/test.fa", "testdata/tmp.noMask.txt", "testdata/expected.noMask.txt", 1, 500, true, 50},
	{"testdata/test.fa", "testdata/tmp.lineLength.txt", "testdata/expected.lineLength.txt", 1, 500, false, 100},
	{"testdata/test.fa", "testdata/tmp.short.txt", "testdata/expected.short.txt", 350, 400, false, 50},
}

func TestMultFaVisualizer(t *testing.T) {
	var err error
	for _, v := range MultFaVisualizerTests {
		multFaVisualizer(v.Infile, v.Outfile, v.Start, v.End, v.NoMask, v.Linelength)
		if !fileio.AreEqual(v.Outfile, v.ExpectedFile) {
			t.Errorf("Error in MultFaVisualizer. Output does not match expected.")
		} else {
			err = os.Remove(v.Outfile)
			exception.PanicOnErr(err)
		}
	}
}
