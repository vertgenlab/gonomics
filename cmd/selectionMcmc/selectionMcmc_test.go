package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"os"
	"testing"
)

var SelectionMcmcTests = []struct {
	VcfFile                 string
	OutFile                 string
	ExpectedFile            string
	Iterations              int
	MuZero                  float64
	SigmaZero               float64
	MuStep                  float64
	SigmaStep               float64
	RandSeed                bool
	SetSeed                 int64
	UnPolarized             bool
	DivergenceAscertainment bool
	FixedSigma              bool
	IntegralError           float64
	Verbose                 int
}{
	{"testdata/N100.S20.AlphaMinus10.Seed20.vcf",
		"testdata/tmp.trace.txt",
		"testdata/expected.trace.txt",
		100,
		-5,
		0.1,
		0.2,
		5,
		false,
		-1,
		false,
		false,
		false,
		1e-7,
		0},
}

func TestSelectionMcmc(t *testing.T) {
	var err error
	var s popgen.McmcSettings
	for _, v := range SelectionMcmcTests {
		s = popgen.McmcSettings{
			Iterations:              v.Iterations,
			MuStep:                  v.MuStep,
			MuZero:                  v.MuZero,
			SigmaStep:               v.SigmaStep,
			SigmaZero:               v.SigmaZero,
			RandSeed:                v.RandSeed,
			SetSeed:                 v.SetSeed,
			UnPolarized:             v.UnPolarized,
			DivergenceAscertainment: v.DivergenceAscertainment,
			FixedSigma:              v.FixedSigma,
			D:                       1,
			IntegralError:           v.IntegralError,
			Verbose:                 v.Verbose,
		}
		selectionMcmc(v.VcfFile, v.OutFile, s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in selectionMcmc. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
