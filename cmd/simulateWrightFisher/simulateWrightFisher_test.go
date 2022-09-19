package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
)

var SimulateWFTests = []struct {
	OutFile       string
	ExpectedFile  string
	PopSize       int
	MutRate       float64
	NumGen        int
	GenomeSize    int
	RFitness      float64
	GcContent     float64
	InitFreq      string
	FitnessString string
	SetSeed       int64
	Verbose       bool
	Fasta         bool
	Vcf           bool
}{
	{"testdata/tmp_without_initFreq.tsv",
		"testdata/expected_without_initFreq.tsv",
		1000,
		1e-4,
		1000,
		1,
		1.02,
		0.5,
		"",
		"",
		5,
		false,
		false,
		false,
	},
	{"testdata/tmp_with_initFreq.tsv",
		"testdata/expected_with_initFreq.tsv",
		1000,
		1e-9,
		500,
		1,
		2,
		0.5,
		"0.25,0.25,0.25,0.25,A",
		"",
		10,
		false,
		false,
		false,
	},
	{"testdata/tmp_with_fitnessString.tsv",
		"testdata/expected_with_fitnessString.tsv",
		1000,
		1e-9,
		500,
		1,
		2,
		0.5,
		"0.25,0.25,0.25,0.25,A",
		"1,1.05,0.95,0.95,A",
		20,
		false,
		false,
		false,
	},
}

func TestSimulateWF(t *testing.T) {
	var err error
	var s popgen.WrightFisherSettings
	for _, v := range SimulateWFTests {
		s = popgen.WrightFisherSettings{
			PopSize:       v.PopSize,
			MutRate:       v.MutRate,
			NumGen:        v.NumGen,
			GenomeSize:    v.GenomeSize,
			RFitness:      v.RFitness,
			GcContent:     v.GcContent,
			InitFreq:      v.InitFreq,
			FitnessString: v.FitnessString,
			SetSeed:       v.SetSeed,
			Verbose:       v.Verbose,
			Fasta:         v.Fasta,
			Vcf:           v.Vcf,
		}
		simulateWrightFisher(v.OutFile, s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in SimulateWrightFisher. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
