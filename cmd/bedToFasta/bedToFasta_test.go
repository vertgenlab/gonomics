package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var BedToFastaTests = []struct {
	FastaFile    string
	BedFile      string
	OutFile      string
	RevComp      bool
	ExpectedFile string
}{
	{"testdata/test.fa", "testdata/test.bed", "testdata/out.fa", false, "testdata/expected.fa"},
	{"testdata/test.fa", "testdata/test.bed", "testdata/outRevComp.fa", true, "testdata/expectedRevComp.fa"},
}

func TestBedToFasta(t *testing.T) {
	var err error
	for _, v := range BedToFastaTests {
		bedToFasta(v.FastaFile, v.BedFile, v.OutFile, v.RevComp)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in bedToFasta. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
