package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
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
	NoHeader             bool
}{
	{InWigFile: "testdata/toTrainSet/toTrainSet.wig",
		InFastaFile:          "testdata/toTrainSet/toTrainSet.fa",
		TrainFile:            "testdata/toTrainSet/tmp.train.txt",
		ValidateFile:         "testdata/toTrainSet/tmp.validate.txt",
		TestFile:             "testdata/toTrainSet/tmp.test.txt",
		WindowSize:           3,
		Stride:               3,
		ValidationProp:       0.1,
		TestingProp:          0.1,
		SetSeed:              5,
		Missing:              -10,
		ExpectedTrainFile:    "testdata/toTrainSet/expected.train.txt",
		ExpectedValidateFile: "testdata/toTrainSet/expected.validate.txt",
		ExpectedTestFile:     "testdata/toTrainSet/expected.test.txt",
		LogTransform:         false,
		IncludeRevComp:       false,
		NoHeader:             false,
	},
	{InWigFile: "testdata/toTrainSet/toTrainSet.wig",
		InFastaFile:          "testdata/toTrainSet/toTrainSet.fa",
		TrainFile:            "testdata/toTrainSet/tmp.log.train.txt",
		ValidateFile:         "testdata/toTrainSet/tmp.log.validate.txt",
		TestFile:             "testdata/toTrainSet/tmp.log.test.txt",
		WindowSize:           3,
		Stride:               3,
		ValidationProp:       0.1,
		TestingProp:          0.1,
		SetSeed:              5,
		Missing:              -10,
		ExpectedTrainFile:    "testdata/toTrainSet/expected.log.train.txt",
		ExpectedValidateFile: "testdata/toTrainSet/expected.log.validate.txt",
		ExpectedTestFile:     "testdata/toTrainSet/expected.log.test.txt",
		LogTransform:         true,
		IncludeRevComp:       false,
		NoHeader:             false,
	},
	{InWigFile: "testdata/toTrainSet/toTrainSet.wig",
		InFastaFile:          "testdata/toTrainSet/toTrainSet.fa",
		TrainFile:            "testdata/toTrainSet/tmp.revComp.train.txt",
		ValidateFile:         "testdata/toTrainSet/tmp.revComp.validate.txt",
		TestFile:             "testdata/toTrainSet/tmp.revComp.test.txt",
		WindowSize:           3,
		Stride:               3,
		ValidationProp:       0.1,
		TestingProp:          0.1,
		SetSeed:              5,
		Missing:              -10,
		ExpectedTrainFile:    "testdata/toTrainSet/expected.revComp.train.txt",
		ExpectedValidateFile: "testdata/toTrainSet/expected.revComp.validate.txt",
		ExpectedTestFile:     "testdata/toTrainSet/expected.revComp.test.txt",
		LogTransform:         false,
		IncludeRevComp:       true,
		NoHeader:             false,
	},
	{InWigFile: "testdata/toTrainSet/toTrainSet.wig",
		InFastaFile:          "testdata/toTrainSet/toTrainSet.fa",
		TrainFile:            "testdata/toTrainSet/tmp.noHeader.train.txt",
		ValidateFile:         "testdata/toTrainSet/tmp.noHeader.validate.txt",
		TestFile:             "testdata/toTrainSet/tmp.noHeader.test.txt",
		WindowSize:           3,
		Stride:               3,
		ValidationProp:       0.1,
		TestingProp:          0.1,
		SetSeed:              5,
		Missing:              -10,
		ExpectedTrainFile:    "testdata/toTrainSet/expected.noHeader.train.txt",
		ExpectedValidateFile: "testdata/toTrainSet/expected.noHeader.validate.txt",
		ExpectedTestFile:     "testdata/toTrainSet/expected.noHeader.test.txt",
		LogTransform:         false,
		IncludeRevComp:       true,
		NoHeader:             true,
	},
}

func TestWigToTrainingSet(t *testing.T) {
	var err error

	for _, v := range WigToTrainingSetTests {
		s := ToTrainingSetSettings{
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
			NoHeader:       v.NoHeader,
		}
		toTrainingSet(s)
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
