package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
)

var PeaksTests = []struct {
	InWig        string
	ChromSizes   string
	OutBed       string
	ExpectedBed  string
	Threshold    float64
	FindMinima   bool
	DefaultValue float64
}{
	{InWig: "testdata/peaks/in_wig_1.wig",
		ChromSizes:   "testdata/peaks/genome.chrom.sizes",
		OutBed:       "testdata/peaks/tmp.bed",
		ExpectedBed:  "testdata/peaks/out_bed_1.bed",
		Threshold:    20,
		FindMinima:   false,
		DefaultValue: 0,
	},
	{InWig: "testdata/peaks/in_wig_1.wig",
		ChromSizes:   "testdata/peaks/genome.chrom.sizes",
		OutBed:       "testdata/peaks/tmp.Minima.bed",
		ExpectedBed:  "testdata/peaks/expected.minima.bed",
		Threshold:    50,
		FindMinima:   true,
		DefaultValue: 100,
	},
}

func TestPeaks(t *testing.T) {
	var err error
	var s PeakSettings
	for _, v := range PeaksTests {
		s = PeakSettings{
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
			t.Errorf("Error in wigTools peaks, created bed and test bed are not equal.")
		} else {
			err = os.Remove(v.OutBed)
			exception.PanicOnErr(err)
		}
	}
}
