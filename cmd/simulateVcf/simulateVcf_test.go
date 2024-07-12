package main

import (
	"testing"

	"github.com/vertgenlab/gonomics/fileio"
)

var SimulateVcfTests = []struct {
	ExpectedFile    string
	OutFile         string
	Alpha           float64
	NumAlleles      int
	NumSites        int
	SetSeed         int64
	BoundAlpha      float64
	BoundBeta       float64
	BoundMultiplier float64
}{
	{"testdata/expected.vcf.gz", "testdata/out.vcf.gz", 4, 100, 100, 11, 0.001, 0.001, 10000},
}

func TestSimulateVcf(t *testing.T) {
	var s Settings
	for _, v := range SimulateVcfTests {
		s = Settings{
			OutFile:         v.OutFile,
			Alpha:           v.Alpha,
			NumSites:        v.NumSites,
			NumAlleles:      v.NumAlleles,
			SetSeed:         v.SetSeed,
			BoundAlpha:      v.BoundAlpha,
			BoundBeta:       v.BoundBeta,
			BoundMultiplier: v.BoundMultiplier,
		}
		simulateVcf(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in simulateVcf.")
		} else {
			fileio.EasyRemove(v.OutFile)
		}
	}
}
