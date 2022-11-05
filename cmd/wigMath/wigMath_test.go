package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var WigMathTests = []struct {
	InFile string
	OutFile string
	ExpectedFile string
	SubtractFile string
	MovingAverageSmoothing int
}{
	{
		InFile: "testdata/in.wig",
		OutFile: "testdata/tmp.wig",
		ExpectedFile: "testdata/expected.subtract.wig",
		SubtractFile: "testdata/subtract.wig",
		MovingAverageSmoothing: 1,
	},
	{
		InFile: "testdata/unsmooth.wig",
		OutFile: "testdata/tmp.smooth.wig",
		ExpectedFile: "testdata/expected.smooth.wig",
		SubtractFile: "",
		MovingAverageSmoothing: 5,
	},
}

func TestWigMath(t *testing.T) {
	var err error
	var s Settings
	for _, v := range WigMathTests {
		s = Settings{
			InFile: v.InFile,
			OutFile: v.OutFile,
			ElementWiseSubtract: v.SubtractFile,
			MovingAverageSmoothing: v.MovingAverageSmoothing,
		}
		wigMath(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in wigMath. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
