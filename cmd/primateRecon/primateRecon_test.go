package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var PrimateReconTests = []struct {
	Infile string
	Outfile string
	ExpectedFile string
	MessyToN bool
}{
	{"testdata/in.fa", "testdata/tmp.fa", "testdata/expected.fa", true},
}

func TestPrimateRecon(t *testing.T) {
	var err error
	for _, v := range PrimateReconTests {
		primateRecon(v.Infile, v.Outfile, v.MessyToN)
		if !fileio.AreEqual(v.Outfile, v.ExpectedFile) {
			t.Errorf("Error in primateRecon. Output did not match expected.")
		} else {
			err = os.Remove(v.Outfile)
			exception.PanicOnErr(err)
		}
	}
}
