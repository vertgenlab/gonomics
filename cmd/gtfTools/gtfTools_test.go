package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var ToBedTests = []struct {
	InFile         string
	OutFile        string
	ExpectedFile   string
	Tss            bool
	ChromSizesFile string
	Merge          bool
}{
	{InFile: "testdata/test.gtf",
		OutFile:        "testdata/tmp.bed",
		ExpectedFile:   "testdata/testOut.bed",
		Tss:            false,
		ChromSizesFile: "",
		Merge:          false,
	},
	{InFile: "testdata/test.gtf",
		OutFile:        "testdata/tmp.tss.bed",
		ExpectedFile:   "testdata/expected.tss.bed",
		Tss:            true,
		ChromSizesFile: "testdata/chr1.chrom.sizes",
		Merge:          false,
	},
}

func TestToBed(t *testing.T) {
	var err error
	var s ToBedSettings
	for _, v := range ToBedTests {
		s = ToBedSettings{
			InFile:        v.InFile,
			OutFile:       v.OutFile,
			Tss:           v.Tss,
			ChromSizeFile: v.ChromSizesFile,
			Merge:         v.Merge,
		}
		toBed(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: in gtfTools toBed, output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
