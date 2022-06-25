package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedFormatTests = []struct {
	InFile         string
	OutFile        string
	ExpectedFile   string
	UCSCToEnsembl  bool
	EnsemblToUCSC  bool
	ScaleNameFloat float64
	PadLength      int
	ChromSizeFile  string
}{
	{"testdata/test.bed", //this test is for scaleNameFloat
		"testdata/test.outFloat.bed",
		"testdata/expected.NameFloat.bed",
		false,
		false,
		10,
		0,
		"",
	},
	{"testdata/test.bed", //this test is for UCSCToEnsembl
		"testdata/test.outEnsembl.bed",
		"testdata/expected.Ensembl.bed",
		true,
		false,
		1,
		0,
		"",
	},
	{"testdata/test.Ensembl.bed", //this test is for UCSCToEnsembl
		"testdata/test.outUCSC.bed",
		"testdata/expected.UCSC.bed",
		false,
		true,
		1,
		0,
		"",
	},
	{"testdata/test.pad.bed",
		"testdata/out.pad.bed",
		"testdata/expected.pad.bed",
		false,
		false,
		1,
		91,
		"testdata/test.chrom.sizes",
	},
}

func TestBedFormat(t *testing.T) {
	var err error
	var s Settings
	for _, v := range BedFormatTests {
		s = Settings{
			InFile:         v.InFile,
			OutFile:        v.OutFile,
			UCSCToEnsembl:  v.UCSCToEnsembl,
			EnsemblToUCSC:  v.EnsemblToUCSC,
			ScaleNameFloat: v.ScaleNameFloat,
			ChromSizeFile:  v.ChromSizeFile,
			PadLength:      v.PadLength,
		}
		bedFormat(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in bedFormat. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
