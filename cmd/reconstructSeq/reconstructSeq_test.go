package main

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

var ReconstructSeqTests = []struct {
	NewickFile           string
	FastaFile            string
	OutFile              string
	ExpectedFile         string
	BiasLeafName         string
	BiasNodeName         string
	NonBiasProbThreshold float64
	BiasN                bool
	HighestProbThreshold float64
	KeepAllSeq           bool
	SubMatrix            bool
	PDnaNode             string
	PDnaNodeMulti        []string
	PDnaOutFile          string
	ExpectedPFile        string
	Precision            float32
}{
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.fa",
		BiasLeafName:         "",
		BiasNodeName:         "",
		NonBiasProbThreshold: 0,
		BiasN:                false,
		HighestProbThreshold: 0,
		KeepAllSeq:           false,
		SubMatrix:            false,
		PDnaNode:             "hca",
		PDnaOutFile:          "testdata/hca1.pfa",
		ExpectedPFile:        "testdata/hca1Expected.pfa",
		Precision:            1e-3,
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.ThresholdPoint8.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.ThresholdPoint8.fa",
		BiasLeafName:         "human",
		BiasNodeName:         "",
		NonBiasProbThreshold: 0.8,
		BiasN:                false,
		HighestProbThreshold: 0,
		KeepAllSeq:           false,
		SubMatrix:            false,
		PDnaNode:             "hga",
		PDnaOutFile:          "testdata/hga1.pfa",
		ExpectedPFile:        "testdata/hga1Expected.pfa",
		Precision:            1e-3,
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.ThresholdPoint99.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.ThresholdPoint99.fa",
		BiasLeafName:         "human",
		BiasNodeName:         "",
		NonBiasProbThreshold: 0.99,
		BiasN:                false,
		HighestProbThreshold: 0,
		KeepAllSeq:           false,
		SubMatrix:            false,
		PDnaNode:             "hoa",
		PDnaOutFile:          "testdata/hoa1.pfa",
		ExpectedPFile:        "testdata/hoa1Expected.pfa",
		Precision:            1e-3,
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.highestProbThreshold99.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.highestProbThreshold99.fa",
		BiasLeafName:         "human",
		BiasNodeName:         "",
		NonBiasProbThreshold: 0,
		BiasN:                false,
		HighestProbThreshold: 0.99,
		KeepAllSeq:           false,
		SubMatrix:            false,
		PDnaNode:             "cba",
		PDnaOutFile:          "testdata/cba1.pfa",
		ExpectedPFile:        "testdata/cba1Expected.pfa",
		Precision:            1e-3,
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.withExtraSeqs.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.withExtraSeqs.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.fa",
		BiasLeafName:         "",
		BiasNodeName:         "",
		NonBiasProbThreshold: 0,
		BiasN:                false,
		HighestProbThreshold: 0,
		KeepAllSeq:           false,
		SubMatrix:            false,
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.withExtraSeqs.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.withExtraSeqs.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.keepAllSeq.fa",
		BiasLeafName:         "",
		BiasNodeName:         "",
		NonBiasProbThreshold: 0,
		BiasN:                false,
		HighestProbThreshold: 0,
		KeepAllSeq:           true,
		SubMatrix:            false,
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.withExtraSeqsRef.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.withExtraSeqsRef.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.keepAllSeqRef.fa",
		BiasLeafName:         "",
		BiasNodeName:         "",
		NonBiasProbThreshold: 0,
		BiasN:                false,
		HighestProbThreshold: 0,
		KeepAllSeq:           true,
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/short.fa",
		OutFile:              "testdata/out.short.biasN.fa",
		ExpectedFile:         "testdata/expected.short.biasN.fa",
		BiasLeafName:         "human",
		BiasNodeName:         "",
		NonBiasProbThreshold: 0.8,
		BiasN:                true,
		HighestProbThreshold: 0,
		KeepAllSeq:           true,
		PDnaNode:             "hga",
		PDnaOutFile:          "testdata/hga2.pfa",
		ExpectedPFile:        "testdata/hga2Expected.pfa",
		Precision:            1e-3,
	},
	{NewickFile: "testdata/allT2T.4d.mod",
		FastaFile:            "testdata/allT2T.fa",
		OutFile:              "testdata/out.allT2T.biasNodeName.fa",
		ExpectedFile:         "testdata/expected.allT2T.biasNodeName.fa",
		BiasLeafName:         "chimpT2Tpri",
		BiasNodeName:         "hcaT2T",
		NonBiasProbThreshold: 0.8,
		BiasN:                false,
		HighestProbThreshold: 0,
		KeepAllSeq:           true,
	},
	{NewickFile: "testdata/4d.genericNames.mod",
		FastaFile:            "testdata/allPossible.oneHuman.fa",
		OutFile:              "testdata/out.AllPossibleOneHuman.fa",
		ExpectedFile:         "testdata/expected.AllPossibleOneHuman.fa",
		BiasLeafName:         "",
		BiasNodeName:         "",
		NonBiasProbThreshold: 0,
		BiasN:                false,
		HighestProbThreshold: 0,
		KeepAllSeq:           false,
		SubMatrix:            false,
		PDnaNodeMulti:        []string{"hca", "hga"},
		PDnaOutFile:          "testdata/multi_hca_hga.pfa",
		ExpectedPFile:        "testdata/multi_hca_hgaExpected.pfa",
		Precision:            1e-3,
	},
}

func TestReconstructSeq(t *testing.T) {
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
			BiasNodeName:         v.BiasNodeName,
			NonBiasProbThreshold: v.NonBiasProbThreshold,
			BiasN:                v.BiasN,
			HighestProbThreshold: v.HighestProbThreshold,
			KeepAllSeq:           v.KeepAllSeq,
			SubMatrix:            v.SubMatrix,
			PDnaNode:             v.PDnaNode,
			PDnaNodeMulti:        v.PDnaNodeMulti,
			PDnaOutFile:          v.PDnaOutFile,
		}
		ReconstructSeq(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: output was not as expected.")
		} else {
			fileio.EasyRemove(v.OutFile)
		}
		if v.PDnaNode != "" || len(v.PDnaNodeMulti) > 0 {
			expectedPFa := pFasta.Read(v.ExpectedPFile)
			reconPFa := pFasta.Read(v.PDnaOutFile)
			if !pFasta.AllAreEqual(expectedPFa, reconPFa, v.Precision) {
				t.Errorf("Error: pfasta output not as expected.")
			} else {
				fileio.EasyRemove(v.PDnaOutFile)
			}
		}
	}
}
