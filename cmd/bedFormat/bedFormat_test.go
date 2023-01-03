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
	LeftPadLength  int
	RightPadLength int
	ChromSizeFile  string
	ToMidpoint     bool
}{
	{InFile: "testdata/test.bed", //this test is for scaleNameFloat
		OutFile:        "testdata/test.outFloat.bed",
		ExpectedFile:   "testdata/expected.NameFloat.bed",
		UCSCToEnsembl:  false,
		EnsemblToUCSC:  false,
		ScaleNameFloat: 10,
		PadLength:      0,
		LeftPadLength:  0,
		RightPadLength: 0,
		ChromSizeFile:  "",
		ToMidpoint:     false,
	},
	{InFile: "testdata/test.bed", //this test is for UCSCToEnsembl
		OutFile:        "testdata/test.outEnsembl.bed",
		ExpectedFile:   "testdata/expected.Ensembl.bed",
		UCSCToEnsembl:  true,
		EnsemblToUCSC:  false,
		ScaleNameFloat: 1,
		PadLength:      0,
		LeftPadLength:  0,
		RightPadLength: 0,
		ChromSizeFile:  "",
		ToMidpoint:     false,
	},
	{InFile: "testdata/test.Ensembl.bed", //this test is for UCSCToEnsembl
		OutFile:        "testdata/test.outUCSC.bed",
		ExpectedFile:   "testdata/expected.UCSC.bed",
		UCSCToEnsembl:  false,
		EnsemblToUCSC:  true,
		ScaleNameFloat: 1,
		PadLength:      0,
		LeftPadLength:  0,
		RightPadLength: 0,
		ChromSizeFile:  "",
		ToMidpoint:     false,
	},
	{InFile: "testdata/pad.bed",
		OutFile:        "testdata/out.pad.bed",
		ExpectedFile:   "testdata/expected.pad.bed",
		UCSCToEnsembl:  false,
		EnsemblToUCSC:  false,
		ScaleNameFloat: 1,
		PadLength:      91,
		LeftPadLength:  0,
		RightPadLength: 0,
		ChromSizeFile:  "testdata/test.chrom.sizes",
		ToMidpoint:     false,
	},
	{InFile: "testdata/test.bed",
		OutFile:        "testdata/test.midpoint.bed",
		ExpectedFile:   "testdata/expected.midpoint.bed",
		UCSCToEnsembl:  false,
		EnsemblToUCSC:  false,
		ScaleNameFloat: 1,
		PadLength:      0,
		LeftPadLength:  0,
		RightPadLength: 0,
		ChromSizeFile:  "testdata/test.chrom.sizes",
		ToMidpoint:     true,
	},
	{InFile: "testdata/pad.bed",
		OutFile:        "testdata/test.leftPad.bed",
		ExpectedFile:   "testdata/expected.leftPad.bed",
		UCSCToEnsembl:  false,
		EnsemblToUCSC:  false,
		ScaleNameFloat: 1,
		PadLength:      0,
		LeftPadLength:  91,
		RightPadLength: 0,
		ChromSizeFile:  "testdata/test.chrom.sizes",
		ToMidpoint:     false,
	},
	{InFile: "testdata/pad.bed",
		OutFile:        "testdata/test.rightPad.bed",
		ExpectedFile:   "testdata/expected.rightPad.bed",
		UCSCToEnsembl:  false,
		EnsemblToUCSC:  false,
		ScaleNameFloat: 1,
		PadLength:      0,
		LeftPadLength:  0,
		RightPadLength: 91,
		ChromSizeFile:  "testdata/test.chrom.sizes",
		ToMidpoint:     false,
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
			EvenPadLength:  v.PadLength,
			LeftPadLength: v.LeftPadLength,
			RightPadLength: v.RightPadLength,
			ToMidpoint:     v.ToMidpoint,
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
