package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var MafFilterTests = []struct {
	Infile       string
	Outfile      string
	Threshold    float64
	ExpectedFile string
}{
	{"testdata/chr22.test.maf", "testdata/tmp.maf", 10000, "testdata/expected.chr22.maf"},
}

func TestMafFilter(t *testing.T) {
	var err error
	for _, v := range MafFilterTests {
		mafFilter(v.Infile, v.Outfile, v.Threshold)
		if !fileio.AreEqual(v.Outfile, v.ExpectedFile) {
			t.Errorf("Error in mafFilter. Output did not match expected.")
		} else {
			err = os.Remove(v.Outfile)
			exception.PanicOnErr(err)
		}
	}
}
