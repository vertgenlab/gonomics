package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var FilterTests = []struct {
	InFile       string
	ChromSizes   string
	OutFile      string
	ExpectedFile string
	Chrom        string
}{
	{InFile: "testdata/in.filter.wig",
		ChromSizes:   "testdata/test.filter.chrom.sizes",
		OutFile:      "testdata/tmp.filter.wig",
		ExpectedFile: "testdata/expected.filter.wig",
		Chrom:        "chr3"},
}

func TestFilter(t *testing.T) {
	var err error
	var s FilterSettings
	for _, v := range FilterTests {
		s = FilterSettings{
			InFile:     v.InFile,
			ChromSizes: v.ChromSizes,
			OutFile:    v.OutFile,
			Chrom:      v.Chrom,
		}
		wigFilter(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in wigTools filter. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
