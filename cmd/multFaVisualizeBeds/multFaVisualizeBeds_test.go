package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var MultFaVisualizeBedsTests = []struct {
	BedFile string
	MultiFaFile string
	OutFormat bool
	NoMask bool
	LineLength int
	OutFiles []string
	ExpectedFiles []string
}{
	{"testdata/test.bed", "testdata/test.fa", false, false, 50, []string{"testdata/chr1_10_100.txt", "testdata/chr1_490_500.txt", "testdata/chr1_5_10.txt"}, []string{"testdata/expected.chr1_10_100.txt", "testdata/expected.chr1_490_500.txt", "testdata/expected.chr1_5_10.txt"}},
	{"testdata/test.bed", "testdata/test.fa", false, false, 100, []string{"testdata/chr1_10_100.txt", "testdata/chr1_490_500.txt", "testdata/chr1_5_10.txt"}, []string{"testdata/expected.long.chr1_10_100.txt", "testdata/expected.chr1_490_500.txt", "testdata/expected.chr1_5_10.txt"}},
}


func TestMultFaVisualizeBeds(t *testing.T) {
	var err error
	for _, v := range MultFaVisualizeBedsTests {
		multFaVisualizeBeds(v.BedFile, v.MultiFaFile, v.OutFormat, v.NoMask, v.LineLength, "testdata/")
		for i := range v.OutFiles {
			if !fileio.AreEqual(v.OutFiles[i], v.ExpectedFiles[i]) {
				t.Errorf("ERror in multFaVisualizeBeds. Output did not match expected.")
			} else {
				err = os.Remove(v.OutFiles[i])
				exception.PanicOnErr(err)
			}
		}
	}
}
