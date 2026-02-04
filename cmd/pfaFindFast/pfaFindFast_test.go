package main

import (
	"math"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var FaFindFastTests = []struct {
	InFile                  string
	OutFile                 string
	ExpectedFile            string
	FirstQueryName          string
	SecondQueryName         string
	WindowSize              int
	RefChromName            string
	RemoveN                 bool
	DivergenceRate          float64
	LongOutput              bool
	OutputAlnPos            bool
	BaseDistToDivThreshold  float64
	BaseDotToSubstThreshold float64
	ConfidentThreshold      float32
}{
	{InFile: "testdata/human_hca_hga.pfa",
		OutFile:                 "testdata/tmp.out.bed",
		ExpectedFile:            "testdata/expected.bed",
		FirstQueryName:          "hca",
		SecondQueryName:         "hga",
		WindowSize:              10,
		RefChromName:            "chr1",
		RemoveN:                 false,
		DivergenceRate:          math.MaxFloat64,
		LongOutput:              false,
		OutputAlnPos:            false,
		BaseDistToDivThreshold:  0.7,
		BaseDotToSubstThreshold: 0.8,
		ConfidentThreshold:      0.8,
	},
	{InFile: "testdata/human_hca_hga.pfa",
		OutFile:                 "testdata/tmp.out.bed",
		ExpectedFile:            "testdata/expected.longOutput.bed",
		FirstQueryName:          "hca",
		SecondQueryName:         "hga",
		WindowSize:              10,
		RefChromName:            "chr1",
		RemoveN:                 false,
		DivergenceRate:          math.MaxFloat64,
		LongOutput:              true,
		OutputAlnPos:            false,
		BaseDistToDivThreshold:  0.7,
		BaseDotToSubstThreshold: 0.8,
		ConfidentThreshold:      0.8,
	},
}

func TestFaFindFast(t *testing.T) {
	var err error
	for _, v := range FaFindFastTests {
		s := Settings{
			InFile:                  v.InFile,
			OutFile:                 v.OutFile,
			FirstQueryName:          v.FirstQueryName,
			SecondQueryName:         v.SecondQueryName,
			WindowSize:              v.WindowSize,
			RefChromName:            v.RefChromName,
			RemoveN:                 v.RemoveN,
			LongOutput:              v.LongOutput,
			DivergenceRate:          v.DivergenceRate,
			OutputAlnPos:            v.OutputAlnPos,
			BaseDistToDivThreshold:  v.BaseDistToDivThreshold,
			BaseDotToSubstThreshold: v.BaseDotToSubstThreshold,
			ConfidentThreshold:      v.ConfidentThreshold,
		}
		pfaFindFast(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in faFindFast. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
