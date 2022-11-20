package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)


var PredictSetToHeatmapTests = []struct {
	Mode string
	InFile       string
	OutFile      string
	ExpectedFile string
	WindowSize   int
	Stride       int
	WithRevComp  bool
	SatMutagenesisBeds string
}{
	{Mode: "PredictSetToHeatmap",
		InFile: "testdata/satMut.predicted.txt",
		OutFile:      "chr1.0.5.txt",
		ExpectedFile: "testdata/expected.heatmap.txt",
		WindowSize:   10,
		Stride:       1,
		WithRevComp:  false,
		SatMutagenesisBeds: "",
	},
}

func TestPredictSetToHeatmap(t *testing.T) {
	var err error
	for _, v := range PredictSetToHeatmapTests {
		s := Settings{
			Mode: v.Mode,
			InFile:      v.InFile,
			OutFile:     v.OutFile,
			WindowSize:  v.WindowSize,
			Stride:      v.Stride,
			WithRevComp: v.WithRevComp,
			SatMutagenesisBeds: v.SatMutagenesisBeds,
		}
		PredictSetToHeatmap(s)
		if !fileio.AreEqual(v.ExpectedFile, v.OutFile) {
			t.Errorf("Error in faToSatMutagenesisPredictSet. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var FaToPredictSetCompleteTests = []struct {
	Mode string
	InFile       string
	OutFile      string
	ExpectedFile string
	WindowSize   int
	Stride       int
	WithRevComp  bool
	SatMutagenesisBeds string
}{
	{Mode: "ToPredictSet",
		InFile: "testdata/test.fa",
		OutFile:      "testdata/tmp.txt",
		ExpectedFile: "testdata/expected.txt",
		WindowSize:   10,
		Stride:       1,
		WithRevComp:  false,
		SatMutagenesisBeds: "",
	},
	{Mode: "ToPredictSet",
		InFile: "testdata/test.fa",
		OutFile:      "testdata/tmp.withRevComp.txt",
		ExpectedFile: "testdata/expected.withRevComp.txt",
		WindowSize:   10,
		Stride:       1,
		WithRevComp:  true,
		SatMutagenesisBeds: "",
	},
}

func TestFaToPredictSetComplete(t *testing.T) {
	var err error
	for _, v := range FaToPredictSetCompleteTests {
		s := Settings{
			Mode: v.Mode,
			InFile:      v.InFile,
			OutFile:     v.OutFile,
			WindowSize:  v.WindowSize,
			Stride:      v.Stride,
			WithRevComp: v.WithRevComp,
			SatMutagenesisBeds: v.SatMutagenesisBeds,
		}
		FaToPredictSetComplete(s)
		if !fileio.AreEqual(v.ExpectedFile, v.OutFile) {
			t.Errorf("Error in faToSatMutagenesisPredictSet. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}


var FaToSatMutagenesisPredictSetTests = []struct {
	Mode string
	InFile       string
	OutFile      string
	ExpectedFile string
	WindowSize   int
	Stride       int
	WithRevComp  bool
	SatMutagenesisBeds string
}{
	{Mode: "ToPredictSet",
		InFile: "testdata/test.fa",
		OutFile:      "testdata/tmp.satMut.txt",
		ExpectedFile: "testdata/expected.satMut.txt",
		WindowSize:   10,
		Stride:       1,
		WithRevComp:  false,
		SatMutagenesisBeds: "testdata/satMut.bed",
	},
}

func TestfaToSatMutagenesisPredictSet(t *testing.T) {
	var err error
	for _, v := range FaToSatMutagenesisPredictSetTests {
		s := Settings{
			Mode: v.Mode,
			InFile:      v.InFile,
			OutFile:     v.OutFile,
			WindowSize:  v.WindowSize,
			Stride:      v.Stride,
			WithRevComp: v.WithRevComp,
			SatMutagenesisBeds: v.SatMutagenesisBeds,
		}
		FaToSatMutagenesisPredictSet(s)
		if !fileio.AreEqual(v.ExpectedFile, v.OutFile) {
			t.Errorf("Error in faToSatMutagenesisPredictSet. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
