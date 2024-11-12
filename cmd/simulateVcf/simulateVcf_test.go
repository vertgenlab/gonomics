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
	SetSeed         int64
	BoundAlpha      float64
	BoundBeta       float64
	BoundMultiplier float64
	RefFile         string
	HasRef          bool
}{
	{
		ExpectedFile:    "testdata/expected.vcf",
		OutFile:         "testdata/out.vcf",
		Alpha:           4,
		NumAlleles:      100,
		NumSites:        100,
		SetSeed:         11,
		BoundAlpha:      0.001,
		BoundBeta:       0.001,
		BoundMultiplier: 10000,
		RefFile:         "",
		HasRef:          false,
	},
	{
		ExpectedFile:    "testdata/expected_2.vcf",
		OutFile:         "testdata/out_2.vcf",
		Alpha:           4,
		NumAlleles:      100,
		NumSites:        10,
		SetSeed:         11,
		BoundAlpha:      0.001,
		BoundBeta:       0.001,
		BoundMultiplier: 10000,
		RefFile:         "testdata/refFa_short.fasta",
		HasRef:          true,
	},
	{
		ExpectedFile:    "testdata/expected_3.vcf",
		OutFile:         "testdata/out_3.vcf",
		Alpha:           4,
		NumAlleles:      100,
		NumSites:        20,
		SetSeed:         29,
		BoundAlpha:      0.001,
		BoundBeta:       0.001,
		BoundMultiplier: 10000,
		RefFile:         "testdata/refFa_short.fasta",
		HasRef:          true,
	},
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
			SetSeed:         v.SetSeed,
			BoundAlpha:      v.BoundAlpha,
			BoundBeta:       v.BoundBeta,
			BoundMultiplier: v.BoundMultiplier,
			RefFile:         v.RefFile,
			HasRef:          v.HasRef,
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
