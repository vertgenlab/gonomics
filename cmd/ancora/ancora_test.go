package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var ancoraBuildTests = []struct {
	SamFileName           string
	RefFile               string
	OutFileA              string
	OutFileB              string
	MultiFaDir            string
	qNameA                string
	qNameB                string
	Delta                 float64
	Gamma                 float64
	Epsilon               float64
	Kappa                 float64
	LikelihoodCacheSize   int
	SetSeed               int64
	ExpectedOutFileA      string
	ExpectedOutFileB      string
	OutMultiFaFileNames   []string
	ExpectedMultiFaFiles  []string
	Verbose               int
	EmpiricalPrior        string
	ProblematicRegionsBed string
}{
	{SamFileName: "testdata/diverged.RefAln.sorted.bam",
		RefFile:              "testdata/ref.fa",
		OutFileA:             "testdata/tmp.OutFileA.fa",
		OutFileB:             "testdata/tmp.OutFileB.fa",
		MultiFaDir:           "testdata/multiFa",
		qNameA:               "Rand_Con_A",
		qNameB:               "Rand_Con_B",
		Delta:                0.01,
		Gamma:                3,
		Epsilon:              0.01,
		Kappa:                0.5,
		LikelihoodCacheSize:  100,
		SetSeed:              19,
		ExpectedOutFileA:     "testdata/expected.OutFileA.fa",
		ExpectedOutFileB:     "testdata/expected.OutFileB.fa",
		OutMultiFaFileNames:  []string{"testdata/multiFa/Sequence_0.fa", "testdata/multiFa/Sequence_1.fa"},
		ExpectedMultiFaFiles: []string{"testdata/multiFa/expected.Sequence_0.fa", "testdata/multiFa/expected.Sequence_1.fa"},
		Verbose:              0,
	},
	{SamFileName: "testdata/diverged.RefAln.sorted.bam",
		RefFile:             "testdata/ref.fa",
		OutFileA:            "testdata/tmp.empirical.OutFileA.fa",
		OutFileB:            "testdata/tmp.empirical.OutFileB.fa",
		qNameA:              "Rand_Con_A",
		qNameB:              "Rand_Con_B",
		Delta:               0.01,
		Epsilon:             0.01,
		Kappa:               0.5,
		LikelihoodCacheSize: 100,
		SetSeed:             19,
		ExpectedOutFileA:    "testdata/expected.empirical.OutFileA.fa",
		ExpectedOutFileB:    "testdata/expected.empirical.OutFileB.fa",
		Verbose:             0,
		EmpiricalPrior:      "testdata/expected.AncoraPrior.txt",
	},
	{SamFileName: "testdata/diverged.RefAln.sorted.bam",
		RefFile:               "testdata/ref.fa",
		OutFileA:              "testdata/tmp.problematic.OutFileA.fa",
		OutFileB:              "testdata/tmp.problematic.OutFileB.fa",
		MultiFaDir:            "testdata/problematicMultiFa",
		qNameA:                "Rand_Con_A",
		qNameB:                "Rand_Con_B",
		Delta:                 0.01,
		Gamma:                 3,
		Epsilon:               0.01,
		Kappa:                 0.5,
		LikelihoodCacheSize:   100,
		SetSeed:               19,
		ExpectedOutFileA:      "testdata/expected.problematic.OutFileA.fa",
		ExpectedOutFileB:      "testdata/expected.problematic.OutFileB.fa",
		OutMultiFaFileNames:   []string{"testdata/problematicMultiFa/Sequence_0.fa", "testdata/problematicMultiFa/Sequence_1.fa"},
		ExpectedMultiFaFiles:  []string{"testdata/problematicMultiFa/expected.Sequence_0.fa", "testdata/problematicMultiFa/expected.Sequence_1.fa"},
		Verbose:               0,
		ProblematicRegionsBed: "testdata/problematicRegion.bed",
	},
}

func TestAncoraBuild(t *testing.T) {
	var err error
	var s BuildSettings
	var i int
	for _, v := range ancoraBuildTests {
		s = BuildSettings{
			SamFileName:           v.SamFileName,
			RefFile:               v.RefFile,
			OutFileA:              v.OutFileA,
			OutFileB:              v.OutFileB,
			MultiFaDir:            v.MultiFaDir,
			qNameA:                v.qNameA,
			qNameB:                v.qNameB,
			Delta:                 v.Delta,
			Gamma:                 v.Gamma,
			Epsilon:               v.Epsilon,
			Kappa:                 v.Kappa,
			LikelihoodCacheSize:   v.LikelihoodCacheSize,
			SetSeed:               v.SetSeed,
			Verbose:               v.Verbose,
			EmpiricalPrior:        v.EmpiricalPrior,
			ProblematicRegionsBed: v.ProblematicRegionsBed,
		}
		ancoraBuild(s)
		if !fileio.AreEqual(v.OutFileA, v.ExpectedOutFileA) {
			t.Errorf("Error in ancora build. Outfile A was not as expected.\n")
		} else {
			err = os.Remove(v.OutFileA)
			exception.PanicOnErr(err)
		}
		if !fileio.AreEqual(v.OutFileB, v.ExpectedOutFileB) {
			t.Errorf("Error in ancora build. Outfile B was not as expected.\n")
		} else {
			err = os.Remove(v.OutFileB)
			exception.PanicOnErr(err)
		}
		for i = range v.OutMultiFaFileNames {
			if !fileio.AreEqual(v.OutMultiFaFileNames[i], v.ExpectedMultiFaFiles[i]) {
				t.Errorf("Error: ancora build. MultiFaFile: %v did not match expected: %v.\n", v.OutMultiFaFileNames[i], v.ExpectedMultiFaFiles[i])
			}
		}
	}
}

var ScoreTests = []struct {
	ScoreType    string
	InFileList   string
	OutFile      string
	ExpectedFile string
}{
	{ScoreType: "baseMatrix",
		InFileList:   "testdata/score/fileList.txt",
		OutFile:      "testdata/score/test.baseMatrix.txt",
		ExpectedFile: "testdata/score/expected.baseMatrix.txt"},
	{ScoreType: "baseMatrixByRefBase",
		InFileList:   "testdata/score/fileList.txt",
		OutFile:      "testdata/score/test.baseMatrixByRefBase.txt",
		ExpectedFile: "testdata/score/expected.baseMatrixByRefBase.txt"},
}

func TestAncoraScore(t *testing.T) {
	var err error
	var s ScoreSettings
	for _, v := range ScoreTests {
		s = ScoreSettings{
			ScoreType:  v.ScoreType,
			InFileList: v.InFileList,
			OutFile:    v.OutFile,
		}
		ancoraScore(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: ancora score. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var PriorTests = []struct {
	SamFileName         string
	ReferenceFile       string
	OutFile             string
	Epsilon             float64
	LikelihoodCacheSize int
	ExpectedFile        string
	PseudoCount         float64
	AsCounts            bool
	MinCoverage         int
}{
	{SamFileName: "testdata/diverged.RefAln.sorted.bam",
		ReferenceFile:       "testdata/ref.fa",
		OutFile:             "testdata/out.AncoraPrior.txt",
		Epsilon:             0.01,
		LikelihoodCacheSize: 100,
		ExpectedFile:        "testdata/expected.AncoraPrior.txt",
		PseudoCount:         0.1,
		AsCounts:            false,
		MinCoverage:         0,
	},
	{SamFileName: "testdata/diverged.RefAln.sorted.bam",
		ReferenceFile:       "testdata/ref.fa",
		OutFile:             "testdata/out.AncoraPrior.AsCounts.txt",
		Epsilon:             0.01,
		LikelihoodCacheSize: 100,
		ExpectedFile:        "testdata/expected.AncoraPrior.AsCounts.txt",
		PseudoCount:         0.1,
		AsCounts:            true,
		MinCoverage:         0,
	},
	{SamFileName: "testdata/diverged.RefAln.sorted.bam",
		ReferenceFile:       "testdata/ref.fa",
		OutFile:             "testdata/out.AncoraPrior.minCoverage.txt",
		Epsilon:             0.01,
		LikelihoodCacheSize: 100,
		ExpectedFile:        "testdata/expected.AncoraPrior.minCoverage.txt",
		PseudoCount:         0.1,
		AsCounts:            false,
		MinCoverage:         30,
	},
}

func TestAncoraPrior(t *testing.T) {
	var s PriorSettings
	var err error
	for _, v := range PriorTests {
		s = PriorSettings{
			SamFileName:         v.SamFileName,
			ReferenceFile:       v.ReferenceFile,
			OutFile:             v.OutFile,
			Epsilon:             v.Epsilon,
			LikelihoodCacheSize: v.LikelihoodCacheSize,
			PseudoCount:         v.PseudoCount,
			AsCounts:            v.AsCounts,
			MinCoverage:         v.MinCoverage,
		}
		ancoraPrior(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in ancora prior. Output was not as expected.\n")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
