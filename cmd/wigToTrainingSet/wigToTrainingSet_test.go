package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var WigToTrainingSetTests = []struct {
	InWigFile            string
	InFastaFile          string
	TrainFile            string
	ValidateFile         string
	TestFile             string
	WindowSize           int
	Stride               int
	ValidationProp       float64
	TestingProp          float64
	SetSeed              int64
	Missing              float64
	ExpectedTrainFile    string
	ExpectedValidateFile string
	ExpectedTestFile     string
	LogTransform         bool
	IncludeRevComp       bool
	NoHeader bool
}{
	{InWigFile: "testdata/in.wig",
		InFastaFile:          "testdata/genome.fa",
		TrainFile:            "testdata/tmp.train.txt",
		ValidateFile:         "testdata/tmp.validate.txt",
		TestFile:             "testdata/tmp.test.txt",
		WindowSize:           3,
		Stride:               3,
		ValidationProp:       0.1,
		TestingProp:          0.1,
		SetSeed:              5,
		Missing:              -10,
		ExpectedTrainFile:    "testdata/expected.train.txt",
		ExpectedValidateFile: "testdata/expected.validate.txt",
		ExpectedTestFile:     "testdata/expected.test.txt",
		LogTransform:         false,
		IncludeRevComp:       false,
		NoHeader: false,
	},
	{InWigFile: "testdata/in.wig",
		InFastaFile:          "testdata/genome.fa",
		TrainFile:            "testdata/tmp.log.train.txt",
		ValidateFile:         "testdata/tmp.log.validate.txt",
		TestFile:             "testdata/tmp.log.test.txt",
		WindowSize:           3,
		Stride:               3,
		ValidationProp:       0.1,
		TestingProp:          0.1,
		SetSeed:              5,
		Missing:              -10,
		ExpectedTrainFile:    "testdata/expected.log.train.txt",
		ExpectedValidateFile: "testdata/expected.log.validate.txt",
		ExpectedTestFile:     "testdata/expected.log.test.txt",
		LogTransform:         true,
		IncludeRevComp:       false,
		NoHeader: false,
	},
	{InWigFile: "testdata/in.wig",
		InFastaFile:          "testdata/genome.fa",
		TrainFile:            "testdata/tmp.revComp.train.txt",
		ValidateFile:         "testdata/tmp.revComp.validate.txt",
		TestFile:             "testdata/tmp.revComp.test.txt",
		WindowSize:           3,
		Stride:               3,
		ValidationProp:       0.1,
		TestingProp:          0.1,
		SetSeed:              5,
		Missing:              -10,
		ExpectedTrainFile:    "testdata/expected.revComp.train.txt",
		ExpectedValidateFile: "testdata/expected.revComp.validate.txt",
		ExpectedTestFile:     "testdata/expected.revComp.test.txt",
		LogTransform:         false,
		IncludeRevComp:       true,
		NoHeader: false,
	},
	{InWigFile: "testdata/in.wig",
		InFastaFile:          "testdata/genome.fa",
		TrainFile:            "testdata/tmp.noHeader.train.txt",
		ValidateFile:         "testdata/tmp.noHeader.validate.txt",
		TestFile:             "testdata/tmp.noHeader.test.txt",
		WindowSize:           3,
		Stride:               3,
		ValidationProp:       0.1,
		TestingProp:          0.1,
		SetSeed:              5,
		Missing:              -10,
		ExpectedTrainFile:    "testdata/expected.noHeader.train.txt",
		ExpectedValidateFile: "testdata/expected.noHeader.validate.txt",
		ExpectedTestFile:     "testdata/expected.noHeader.test.txt",
		LogTransform:         false,
		IncludeRevComp:       true,
		NoHeader: true,
	},
}

func TestWigToTrainingSet(t *testing.T) {
	var err error

	for _, v := range WigToTrainingSetTests {
		s := Settings{
			InWigFile:      v.InWigFile,
			InFastaFile:    v.InFastaFile,
			TrainFile:      v.TrainFile,
			ValidateFile:   v.ValidateFile,
			TestFile:       v.TestFile,
			WindowSize:     v.WindowSize,
			Stride:         v.Stride,
			ValidationProp: v.ValidationProp,
			TestingProp:    v.TestingProp,
			SetSeed:        v.SetSeed,
			Missing:        v.Missing,
			LogTransform:   v.LogTransform,
			IncludeRevComp: v.IncludeRevComp,
			NoHeader: v.NoHeader,
		}
		wigToTrainingSet(s)
		if !fileio.AreEqual(v.TrainFile, v.ExpectedTrainFile) {
			t.Errorf("Error in WigToTrainingSet. Training file was not as expected.")
		} else {
			err = os.Remove(v.TrainFile)
			exception.PanicOnErr(err)
		}

		if !fileio.AreEqual(v.TestFile, v.ExpectedTestFile) {
			t.Errorf("Error in WigToTrainingSet. Test file was not as expected.")
		} else {
			err = os.Remove(v.TestFile)
			exception.PanicOnErr(err)
		}

		if !fileio.AreEqual(v.ValidateFile, v.ExpectedValidateFile) {
			t.Errorf("Error in WigToTrainingSet. Validation file was not as expected.")
		} else {
			err = os.Remove(v.ValidateFile)
			exception.PanicOnErr(err)
		}
	}
}
