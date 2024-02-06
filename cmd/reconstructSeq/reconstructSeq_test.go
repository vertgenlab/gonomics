package main

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
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
	KeepAllSeq           bool
	OutPfaFile           string
	PfaNamesSplit        []string // directly input []string. Bypass the splitting of 1 string
	ExpectedPfaFile      string
}{
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.fa",
		BiasLeafName:         "",
		NonBiasProbThreshold: 0,
		HighestProbThreshold: 0,
		KeepAllSeq:           false,
		OutPfaFile:           "",
		PfaNamesSplit:        []string{""},
		ExpectedPfaFile:      "",
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.ThresholdPoint8.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.ThresholdPoint8.fa",
		BiasLeafName:         "human",
		NonBiasProbThreshold: 0.8,
		HighestProbThreshold: 0,
		KeepAllSeq:           false,
		OutPfaFile:           "",
		PfaNamesSplit:        []string{""},
		ExpectedPfaFile:      "",
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.ThresholdPoint99.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.ThresholdPoint99.fa",
		BiasLeafName:         "human",
		NonBiasProbThreshold: 0.99,
		HighestProbThreshold: 0,
		KeepAllSeq:           false,
		OutPfaFile:           "",
		PfaNamesSplit:        []string{""},
		ExpectedPfaFile:      "",
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.highestProbThreshold99.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.highestProbThreshold99.fa",
		BiasLeafName:         "human",
		NonBiasProbThreshold: 0,
		HighestProbThreshold: 0.99,
		KeepAllSeq:           false,
		OutPfaFile:           "",
		PfaNamesSplit:        []string{""},
		ExpectedPfaFile:      "",
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.withExtraSeqs.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.withExtraSeqs.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.fa",
		BiasLeafName:         "",
		NonBiasProbThreshold: 0,
		HighestProbThreshold: 0,
		KeepAllSeq:           false,
		OutPfaFile:           "",
		PfaNamesSplit:        []string{""},
		ExpectedPfaFile:      "",
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.withExtraSeqs.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.withExtraSeqs.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.keepAllSeq.fa",
		BiasLeafName:         "",
		NonBiasProbThreshold: 0,
		HighestProbThreshold: 0,
		KeepAllSeq:           true,
		OutPfaFile:           "",
		PfaNamesSplit:        []string{""},
		ExpectedPfaFile:      "",
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.withExtraSeqsRef.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.withExtraSeqsRef.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.keepAllSeqRef.fa",
		BiasLeafName:         "",
		NonBiasProbThreshold: 0,
		HighestProbThreshold: 0,
		KeepAllSeq:           true,
		OutPfaFile:           "",
		PfaNamesSplit:        []string{""},
		ExpectedPfaFile:      "",
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/short.fa",
		OutFile:              "testdata/out.short.fa",
		ExpectedFile:         "testdata/expected.short.fa",
		BiasLeafName:         "",
		NonBiasProbThreshold: 0,
		HighestProbThreshold: 0,
		KeepAllSeq:           true,
		OutPfaFile:           "testdata/out.short.pFa",
		PfaNamesSplit:        []string{"hca", "hoa"},
		ExpectedPfaFile:      "testdata/expected.short.pFa",
	},
}

func TestReconstructSeq(t *testing.T) {
	var err error
	var hum fasta.Fasta = fasta.Fasta{Name: "human"}
	var chi fasta.Fasta = fasta.Fasta{Name: "chimp"}
	var bon fasta.Fasta = fasta.Fasta{Name: "bonobo"}
	var gor fasta.Fasta = fasta.Fasta{Name: "gorilla"}
	var ora fasta.Fasta = fasta.Fasta{Name: "orangutan"}
	var species []fasta.Fasta
	var h, c, b, g, o dna.Base // the values of the human, humanAlt, chimp, ... bases
	var hBases []dna.Base = []dna.Base{dna.A, dna.N, dna.Gap}
	var possibleBases []dna.Base = []dna.Base{dna.A, dna.C, dna.G, dna.T, dna.N, dna.Gap}

	// given that human is A or N or Gap, we will now go through all possible combinations
	for _, h = range hBases {
		for _, c = range possibleBases {
			for _, b = range possibleBases {
				for _, g = range possibleBases {
					for _, o = range possibleBases {
						hum.Seq = append(hum.Seq, h)
						chi.Seq = append(chi.Seq, c)
						bon.Seq = append(bon.Seq, b)
						gor.Seq = append(gor.Seq, g)
						ora.Seq = append(ora.Seq, o)
					}
				}
			}
		}
	}

	species = []fasta.Fasta{hum, chi, bon, gor, ora}
	fasta.Write("testdata/allPossible.oneHuman.fa", species)

	var s Settings
	for _, v := range ReconstructSeqTests {
		s = Settings{
			NewickInput:          v.NewickFile,
			FastaInput:           v.FastaFile,
			OutFile:              v.OutFile,
			BiasLeafName:         v.BiasLeafName,
			NonBiasProbThreshold: v.NonBiasProbThreshold,
			HighestProbThreshold: v.HighestProbThreshold,
			KeepAllSeq:           v.KeepAllSeq,
			OutPfaFile:           v.OutPfaFile,
			PfaNames:             v.PfaNamesSplit,
		}
		ReconstructSeq(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
		if v.OutPfaFile != "" {
			outPfa := pFasta.Read(v.OutPfaFile)
			expectedpFa := pFasta.Read(v.ExpectedPfaFile)
			// TODO: is this precision level ok?
			if !pFasta.AllAreEqual(outPfa, expectedpFa, 1e-040) { // test fails when precision is set to smallest positive float32 5.4e-079
				t.Errorf("Error: output Pfa was not as expected.")
			} else {
				err = os.Remove(v.OutPfaFile)
				exception.PanicOnErr(err)
			}
		}
	}
}
