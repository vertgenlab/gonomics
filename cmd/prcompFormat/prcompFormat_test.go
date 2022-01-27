package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var PrcompFormatTests = []struct {
	Infile       string
	Outfile      string
	ExpectedFile string
}{
	{"testdata/test.fa", "testdata/tmp.tsv", "testdata/expected.tsv"},
}

func TestPrcompFormat(t *testing.T) {
	var err error
	for _, v := range PrcompFormatTests {
		prcompFormat(v.Infile, v.Outfile, false)
		if !fileio.AreEqual(v.Outfile, v.ExpectedFile) {
			t.Errorf("Error in prcompFormat. Output did not match expected.")
		} else {
			err = os.Remove(v.Outfile)
			exception.PanicOnErr(err)
		}
	}
}
