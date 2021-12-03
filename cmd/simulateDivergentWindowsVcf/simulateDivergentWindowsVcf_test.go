package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var SimulateDivergentWindowsVcfTests = []struct {
	ExpectedUpper string
	ExpectedLower string
	TmpUpper string
	TmpLower string
	Alpha float64
	NumAlleles int
	NumTotalSites int
	NumWindowSites int
	NumWindows int
	RandSeed bool
	SetSeed int64
	BoundAlpha float64
	BoundBeta float64
	BoundMultiplier float64
	UpperPercentile float64
	LowerPercentile float64
}{
	{"testdata/upper.vcf",
		"testdata/lower.vcf",
		"testdata/tmpUpper.vcf",
		"testdata/tmpLower.vcf",
		0.01,
		100,
		1000,
		10,
		100,
		false,
		11,
		0.001,
		0.001,
		10000,
		0.9,
		0.1,
		},
}

func TestSimulateDivergentWindowsVcf(t *testing.T) {
	var err error
	var s Settings
	for _, v := range SimulateDivergentWindowsVcfTests {
		s = Settings {
			UpperOut: v.TmpUpper,
			LowerOut: v.TmpLower,
			Alpha: v.Alpha,
			NumAlleles: v.NumAlleles,
			NumTotalSites: v.NumTotalSites,
			NumWindowSites: v.NumWindowSites,
			NumWindows: v.NumWindows,
			RandSeed: v.RandSeed,
			SetSeed: v.SetSeed,
			BoundAlpha: v.BoundAlpha,
			BoundBeta: v.BoundBeta,
			BoundMultiplier: v.BoundMultiplier,
			UpperPercentile: v.UpperPercentile,
			LowerPercentile: v.LowerPercentile,
		}
		simulateDivergentWindowsVcf(s)
		if !fileio.AreEqual(v.ExpectedUpper, v.TmpUpper) {
			t.Errorf("Error in simulateDivergentWindowsVcf. Upper file does not match expected.")
		} else {
			err = os.Remove(v.TmpUpper)
			exception.PanicOnErr(err)
		}
		if !fileio.AreEqual(v.ExpectedLower, v.TmpLower) {
			t.Errorf("Error in simulateDivergentWindowsVcf. Lower file does not match expected.")
		} else {
			err = os.Remove(v.TmpLower)
			exception.PanicOnErr(err)
		}
	}
}