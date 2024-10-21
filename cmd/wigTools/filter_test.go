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
	{InFile: "testdata/filter/in.filter.wig",
		ChromSizes:   "testdata/filter/test.filter.chrom.sizes",
		OutFile:      "testdata/filter/tmp.filter.wig",
		ExpectedFile: "testdata/filter/expected.filter.wig",
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
