package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var ReconstructSeqTests = []struct {
	NewickFile           string
	FastaFile            string
	OutFile              string
	ExpectedFile         string
	BiasLeafName         string
	NonBiasProbThreshold float64
	HighestProbThreshold float64
}{
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.fa",
		BiasLeafName:         "",
		NonBiasProbThreshold: 0,
		HighestProbThreshold: 0,
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.ThresholdPoint8.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.ThresholdPoint8.fa",
		BiasLeafName:         "human",
		NonBiasProbThreshold: 0.8,
		HighestProbThreshold: 0,
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.ThresholdPoint99.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.ThresholdPoint99.fa",
		BiasLeafName:         "human",
		NonBiasProbThreshold: 0.99,
		HighestProbThreshold: 0,
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.highestProbThreshold99.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.highestProbThreshold99.fa",
		BiasLeafName:         "human",
		NonBiasProbThreshold: 0,
		HighestProbThreshold: 0.99,
	},
}

func TestReconstructSeq(t *testing.T) {
	var err error
	var s Settings
	for _, v := range ReconstructSeqTests {
		s = Settings{
			NewickInput:          v.NewickFile,
			FastaInput:           v.FastaFile,
			OutFile:              v.OutFile,
			BiasLeafName:         v.BiasLeafName,
			NonBiasProbThreshold: v.NonBiasProbThreshold,
			HighestProbThreshold: v.HighestProbThreshold,
		}
		ReconstructSeq(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
