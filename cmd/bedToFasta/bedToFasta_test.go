package main

import (
	"testing"
	"os"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var BedToFastaTests = []struct {
	FastaFile string
	BedFile string
	OutFile string
	ExpectedFile string
}{
	{"testdata/test.fa", "testdata/test.bed", "testdata/out.fa", "testdata/expected.fa"},
}

func TestBedToFasta(t *testing.T) {
	var err error
	for _, v := range BedToFastaTests {
		bedToFasta(v.FastaFile, v.BedFile, v.OutFile)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in bedToFasta. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}