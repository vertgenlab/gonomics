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
	ScalarMultiply float64
	ScalarDivide float64
	AddFile string
	SubtractFile           string
	MovingAverageSmoothing int
	AbsoluteError string
	AbsolutePercentError string
	Missing float64
	Pearson string
	SamplingFrequency float64
}{
	{
		InFile:                 "testdata/in.wig",
		OutFile:                "testdata/tmp.add.wig",
		ExpectedFile:           "testdata/expected.add.wig",
		ScalarMultiply: 1,
		ScalarDivide: 1,
		AddFile: "testdata/second.wig",
		SubtractFile:           "",
		MovingAverageSmoothing: 1,
		AbsoluteError: "",
		AbsolutePercentError: "",
		Missing: -10,
		Pearson: "",
		SamplingFrequency: 1,
	},
	{
		InFile:                 "testdata/in.wig",
		OutFile:                "testdata/tmp.wig",
		ExpectedFile:           "testdata/expected.subtract.wig",
		ScalarMultiply: 1,
		ScalarDivide: 1,
		AddFile: "",
		SubtractFile:           "testdata/second.wig",
		MovingAverageSmoothing: 1,
		AbsoluteError: "",
		AbsolutePercentError: "",
		Missing: -10,
		Pearson: "",
		SamplingFrequency: 1,
	},
	{
		InFile:                 "testdata/unsmooth.wig",
		OutFile:                "testdata/tmp.smooth.wig",
		ExpectedFile:           "testdata/expected.smooth.wig",
		ScalarMultiply: 1,
		ScalarDivide: 1,
		AddFile: "",
		SubtractFile:           "",
		MovingAverageSmoothing: 5,
		AbsoluteError: "",
		AbsolutePercentError: "",
		Missing: -10,
		Pearson: "",
		SamplingFrequency: 1,
	},
	{
		InFile:                 "testdata/in.wig",
		OutFile:                "testdata/tmp.absError.wig",
		ExpectedFile:           "testdata/expected.absError.wig",
		ScalarMultiply: 1,
		ScalarDivide: 1,
		AddFile: "",
		SubtractFile:           "",
		MovingAverageSmoothing: 1,
		AbsoluteError: "testdata/second.wig",
		AbsolutePercentError: "",
		Missing: -10,
		Pearson: "",
		SamplingFrequency: 1,
	},
	{
		InFile:                 "testdata/in.wig",
		OutFile:                "testdata/tmp.absPercentError.wig",
		ExpectedFile:           "testdata/expected.absPercentError.wig",
		ScalarMultiply: 1,
		ScalarDivide: 1,
		AddFile: "",
		SubtractFile:           "",
		MovingAverageSmoothing: 1,
		AbsoluteError: "",
		AbsolutePercentError: "testdata/second.wig",
		Missing: -10,
		Pearson: "",
		SamplingFrequency: 1,
	},
	{
		InFile: "testdata/in.wig",
		OutFile: "testdata/tmp.Pearson.txt",
		ExpectedFile: "testdata/expected.Pearson.txt",
		ScalarMultiply: 1,
		ScalarDivide: 1,
		AddFile: "",
		SubtractFile: "",
		MovingAverageSmoothing: 1,
		AbsoluteError: "",
		AbsolutePercentError: "",
		Missing: -10,
		Pearson: "testdata/second.wig",
		SamplingFrequency: 1,
	},
	{
		InFile: "testdata/in.wig",
		OutFile: "testdata/tmp.Mult50.wig",
		ExpectedFile: "testdata/expected.Mult50.wig",
		ScalarMultiply: 50,
		ScalarDivide: 1,
		AddFile: "",
		SubtractFile: "",
		MovingAverageSmoothing: 1,
		AbsoluteError: "",
		AbsolutePercentError: "",
		Missing: -10,
		Pearson: "",
		SamplingFrequency: 1,
	},
	{
		InFile: "testdata/in.wig",
		OutFile: "testdata/tmp.divide4.wig",
		ExpectedFile: "testdata/expected.divide4.wig",
		ScalarMultiply: 1,
		ScalarDivide: 4,
		AddFile: "",
		SubtractFile: "",
		MovingAverageSmoothing: 1,
		AbsoluteError: "",
		AbsolutePercentError: "",
		Missing: -10,
		Pearson: "",
		SamplingFrequency: 1,
	},
}

func TestWigMath(t *testing.T) {
	var err error
	var s Settings
	for _, v := range WigMathTests {
		s = Settings{
			InFile:                 v.InFile,
			OutFile:                v.OutFile,
			ScalarMultiply: v.ScalarMultiply,
			ScalarDivide: v.ScalarDivide,
			ElementWiseAdd: v.AddFile,
			ElementWiseSubtract:    v.SubtractFile,
			MovingAverageSmoothing: v.MovingAverageSmoothing,
			AbsoluteError: v.AbsoluteError,
			AbsolutePercentError: v.AbsolutePercentError,
			Missing: v.Missing,
			Pearson: v.Pearson,
			SamplingFrequency: v.SamplingFrequency,
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
