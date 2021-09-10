package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var SimulateVcfTests = []struct {
	ExpectedFile    string
	OutFile         string
	Alpha           float64
	NumAlleles      int
	NumSites        int
	RandSeed        bool
	SetSeed         int64
	BoundAlpha      float64
	BoundBeta       float64
	BoundMultiplier float64
}{
	{"testdata/expected.vcf", "testdata/out.vcf", 4, 100, 100, false, 11, 0.001, 0.001, 10000},
}

func TestSimulateVcf(t *testing.T) {
	var err error
	var s Settings
	for _, v := range SimulateVcfTests {
		s = Settings{
			OutFile:         v.OutFile,
			Alpha:           v.Alpha,
			NumSites:        v.NumSites,
			NumAlleles:      v.NumAlleles,
			RandSeed:        v.RandSeed,
			SetSeed:         v.SetSeed,
			BoundAlpha:      v.BoundAlpha,
			BoundBeta:       v.BoundBeta,
			BoundMultiplier: v.BoundMultiplier,
		}
		simulateVcf(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in simulateVcf.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
