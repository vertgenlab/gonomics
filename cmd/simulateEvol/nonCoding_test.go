package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var NonCodingTests = []struct {
	TreeFile               string
	FastaFile              string
	OutFile                string
	SetSeed                int64
	NumNodes               int
	GammaAlpha             float64
	GammaBeta              float64
	GcContent              float64
	LenSeq                 int
	SubstitutionMatrixFile string
	UnitBranchLength       float64
	ExpectedFile           string
	NewickOut              string
	ExpectedNewick         string
}{
	{TreeFile: "",
		FastaFile:              "",
		OutFile:                "testdata/test.NonCoding.fa",
		SetSeed:                29,
		NumNodes:               17,
		GammaAlpha:             1,
		GammaBeta:              50,
		GcContent:              0.41,
		LenSeq:                 50,
		SubstitutionMatrixFile: "",
		UnitBranchLength:       -100,
		NewickOut:              "testdata/test.NewickOut.nh",
		ExpectedFile:           "testdata/expected.NonCoding.fa",
		ExpectedNewick:         "testdata/expected.NewickOut.nh",
	},
	{TreeFile: "testdata/newickLongBranches.txt",
		FastaFile:              "testdata/rand.fa",
		OutFile:                "testdata/test.preMade.fa",
		SetSeed:                31,
		NumNodes:               0,
		GammaAlpha:             1,
		GammaBeta:              50,
		GcContent:              0.41,
		LenSeq:                 50,
		SubstitutionMatrixFile: "testdata/substitutionMatrix.txt",
		UnitBranchLength:       0.5,
		NewickOut:              "testdata/test.NewickOut.PreMade.nh",
		ExpectedNewick:         "testdata/expected.NewickOut.PreMade.nh",
		ExpectedFile:           "testdata/expected.NonCoding.PreMade.fa",
	},
}

func TestNonCoding(t *testing.T) {
	var s NonCodingSettings
	var err error
	for _, v := range NonCodingTests {
		s = NonCodingSettings{
			TreeFile:               v.TreeFile,
			FastaFile:              v.FastaFile,
			OutFile:                v.OutFile,
			UnitBranchLength:       v.UnitBranchLength,
			SubstitutionMatrixFile: v.SubstitutionMatrixFile,
			NumNodes:               v.NumNodes,
			GammaAlpha:             v.GammaAlpha,
			GammaBeta:              v.GammaBeta,
			GcContent:              v.GcContent,
			LenSeq:                 v.LenSeq,
			SetSeed:                v.SetSeed,
			NewickOut:              v.NewickOut,
		}
		NonCoding(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: simulateEvol nonCoding output was not as expected.\n")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
		if !fileio.AreEqual(v.NewickOut, v.ExpectedNewick) {
			t.Errorf("Error: simulateEvol nonCoding output Newick was not as expected.\n")
		} else {
			err = os.Remove(v.NewickOut)
			exception.PanicOnErr(err)
		}
	}
}
