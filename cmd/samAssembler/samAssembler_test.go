package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var SamAssemblerTests = []struct {
	SamFileName          string
	RefFile              string
	OutFileA             string
	OutFileB             string
	MultiFaDir           string
	qNameA               string
	qNameB               string
	Delta                float64
	Gamma                float64
	Epsilon              float64
	Kappa                float64
	LikelihoodCacheSize  int
	SetSeed              int64
	ExpectedOutFileA     string
	ExpectedOutFileB     string
	OutMultiFaFileNames  []string
	ExpectedMultiFaFiles []string
	Verbose              int
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
}

func TestSamAssembler(t *testing.T) {
	var err error
	var s Settings
	var i int
	for _, v := range SamAssemblerTests {
		s = Settings{
			SamFileName:         v.SamFileName,
			RefFile:             v.RefFile,
			OutFileA:            v.OutFileA,
			OutFileB:            v.OutFileB,
			MultiFaDir:          v.MultiFaDir,
			qNameA:              v.qNameA,
			qNameB:              v.qNameB,
			Delta:               v.Delta,
			Gamma:               v.Gamma,
			Epsilon:             v.Epsilon,
			Kappa:               v.Kappa,
			LikelihoodCacheSize: v.LikelihoodCacheSize,
			SetSeed:             v.SetSeed,
			Verbose:             v.Verbose,
		}
		samAssembler(s)
		if !fileio.AreEqual(v.OutFileA, v.ExpectedOutFileA) {
			t.Errorf("Error in samAssembler. Outfile A was not as expected.")
		} else {
			err = os.Remove(v.OutFileA)
			exception.PanicOnErr(err)
		}
		if !fileio.AreEqual(v.OutFileB, v.ExpectedOutFileB) {
			t.Errorf("Error in samAssembler. Outfile B was not as expected.")
		} else {
			err = os.Remove(v.OutFileB)
			exception.PanicOnErr(err)
		}
		for i = range v.OutMultiFaFileNames {
			if !fileio.AreEqual(v.OutMultiFaFileNames[i], v.ExpectedMultiFaFiles[i]) {
				t.Errorf("Error in samAssembler. MultiFaFile: %v did not match expected: %v.\n", v.OutMultiFaFileNames[i], v.ExpectedMultiFaFiles[i])
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

func TestScore(t *testing.T) {
	var err error
	for _, v := range ScoreTests {
		samAssemblerScore(v.ScoreType, v.InFileList, v.OutFile)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: samAssemblerScore. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
