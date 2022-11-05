package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var WigMathTests = []struct {
	InFile                 string
	OutFile                string
	ExpectedFile           string
	SubtractFile           string
	MovingAverageSmoothing int
	AbsoluteError          string
	AbsolutePercentError   string
}{
	{
		InFile:                 "testdata/in.wig",
		OutFile:                "testdata/tmp.wig",
		ExpectedFile:           "testdata/expected.subtract.wig",
		SubtractFile:           "testdata/subtract.wig",
		MovingAverageSmoothing: 1,
		AbsoluteError:          "",
		AbsolutePercentError:   "",
	},
	{
		InFile:                 "testdata/unsmooth.wig",
		OutFile:                "testdata/tmp.smooth.wig",
		ExpectedFile:           "testdata/expected.smooth.wig",
		SubtractFile:           "",
		MovingAverageSmoothing: 5,
		AbsoluteError:          "",
		AbsolutePercentError:   "",
	},
	{
		InFile:                 "testdata/in.wig",
		OutFile:                "testdata/tmp.absError.wig",
		ExpectedFile:           "testdata/expected.absError.wig",
		SubtractFile:           "",
		MovingAverageSmoothing: 1,
		AbsoluteError:          "testdata/subtract.wig",
		AbsolutePercentError:   "",
	},
	{
		InFile:                 "testdata/in.wig",
		OutFile:                "testdata/tmp.absPercentError.wig",
		ExpectedFile:           "testdata/expected.absPercentError.wig",
		SubtractFile:           "",
		MovingAverageSmoothing: 1,
		AbsoluteError:          "",
		AbsolutePercentError:   "testdata/subtract.wig",
	},
}

func TestWigMath(t *testing.T) {
	var err error
	var s Settings
	for _, v := range WigMathTests {
		s = Settings{
			InFile:                 v.InFile,
			OutFile:                v.OutFile,
			ElementWiseSubtract:    v.SubtractFile,
			MovingAverageSmoothing: v.MovingAverageSmoothing,
			AbsoluteError:          v.AbsoluteError,
			AbsolutePercentError:   v.AbsolutePercentError,
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
