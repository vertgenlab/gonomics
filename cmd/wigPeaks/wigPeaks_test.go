package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
)

// use struct to specify testdata.
var wigPeaksTests = []struct {
	InWig        string
	ChromSizes   string
	OutBed       string
	ExpectedBed  string
	Threshold    float64
	FindMinima   bool
	DefaultValue float64
}{
	{InWig: "testdata/in_wig_1.wig",
		ChromSizes:   "testdata/genome.chrom.sizes",
		OutBed:       "testdata/tmp.bed",
		ExpectedBed:  "testdata/out_bed_1.bed",
		Threshold:    20,
		FindMinima:   false,
		DefaultValue: 0,
	},
	{InWig: "testdata/in_wig_1.wig",
		ChromSizes:   "testdata/genome.chrom.sizes",
		OutBed:       "testdata/tmp.Minima.bed",
		ExpectedBed:  "testdata/expected.minima.bed",
		Threshold:    50,
		FindMinima:   true,
		DefaultValue: 100,
	},
}

func TestWigPeaks(t *testing.T) {
	var err error
	var s Settings
	for _, v := range wigPeaksTests {
		s = Settings{
			InWig:        v.InWig,
			ChromSizes:   v.ChromSizes,
			OutBed:       v.OutBed,
			Threshold:    v.Threshold,
			FindMinima:   v.FindMinima,
			DefaultValue: v.DefaultValue,
		}
		wigPeaks(s)
		records := bed.Read(v.OutBed)
		expected := bed.Read(v.ExpectedBed)
		if !bed.AllAreEqual(records, expected) {
			t.Errorf("Error in wigPeaks, created bed and test bed are not equal.")
		} else {
			err = os.Remove(v.OutBed)
			exception.PanicOnErr(err)
		}
	}
}
