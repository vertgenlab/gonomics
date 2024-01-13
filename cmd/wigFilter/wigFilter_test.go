package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var WigFilterTests = []struct {
	InFile       string
	ChromSizes   string
	OutFile      string
	ExpectedFile string
	Chrom        string
}{
	{InFile: "testdata/in.wig",
		ChromSizes:   "testdata/test.chrom.sizes",
		OutFile:      "testdata/tmp.wig",
		ExpectedFile: "testdata/expected.wig",
		Chrom:        "chr3"},
}

func TestWigFilter(t *testing.T) {
	var err error
	var s Settings
	for _, v := range WigFilterTests {
		s = Settings{
			InFile:     v.InFile,
			ChromSizes: v.ChromSizes,
			OutFile:    v.OutFile,
			Chrom:      v.Chrom,
		}
		wigFilter(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in wigFilter. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
