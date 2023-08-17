package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"math"
	"os"
	"testing"
)

var BedpeFilterTests = []struct {
	InFile         string
	OutFile        string
	ExpectedFile   string
	MinScore       int
	MaxScore       int
	MinDistance    int
	MaxDistance    int
	MinStart       int
	MaxStart       int
	Chrom          string
	OnlyInterChrom bool
	OnlyIntraChrom bool
}{
	{InFile: "testdata/testScoreFilter.bedpe",
		OutFile:        "testdata/tmp.bedpe",
		ExpectedFile:   "testdata/expectedScore.bedpe",
		MinScore:       6,
		MaxScore:       math.MaxInt,
		MinDistance:    0,
		MaxDistance:    math.MaxInt,
		MinStart:       0,
		MaxStart:       math.MaxInt,
		Chrom:          "chr1",
		OnlyIntraChrom: false,
		OnlyInterChrom: false},
	{InFile: "testdata/testDistanceFilter.bedpe",
		OutFile:        "testdata/tmpDist.bedpe",
		ExpectedFile:   "testdata/expectedDistance.bedpe",
		MinScore:       6,
		MaxScore:       50,
		MinDistance:    0,
		MaxDistance:    100,
		MinStart:       0,
		MaxStart:       50,
		Chrom:          "chr1",
		OnlyIntraChrom: false,
		OnlyInterChrom: false},
	{InFile: "testdata/testStartFilter.bedpe",
		OutFile:        "testdata/tmpStart.bedpe",
		ExpectedFile:   "testdata/expectedStart.bedpe",
		MinScore:       0,
		MaxScore:       50,
		MinDistance:    0,
		MaxDistance:    100,
		MinStart:       0,
		MaxStart:       30,
		Chrom:          "chr1",
		OnlyIntraChrom: false,
		OnlyInterChrom: false},
	{InFile: "testdata/testChromFilter.bedpe",
		OutFile:        "testdata/tmpChrom.bedpe",
		ExpectedFile:   "testdata/expectedChrom.bedpe",
		MinScore:       0,
		MaxScore:       50,
		MinDistance:    0,
		MaxDistance:    100,
		MinStart:       0,
		MaxStart:       math.MaxInt,
		Chrom:          "chr1",
		OnlyIntraChrom: false,
		OnlyInterChrom: false},
	{InFile: "testdata/testChromFilter.bedpe",
		OutFile:        "testdata/tmpInter.bedpe",
		ExpectedFile:   "testdata/expectedInter.bedpe",
		MinScore:       0,
		MaxScore:       50,
		MinDistance:    0,
		MaxDistance:    100,
		MinStart:       0,
		MaxStart:       math.MaxInt,
		Chrom:          "chr1",
		OnlyIntraChrom: false,
		OnlyInterChrom: true},
	{InFile: "testdata/testChromFilter.bedpe",
		OutFile:        "testdata/tmpIntra.bedpe",
		ExpectedFile:   "testdata/expectedIntra.bedpe",
		MinScore:       0,
		MaxScore:       50,
		MinDistance:    0,
		MaxDistance:    100,
		MinStart:       0,
		MaxStart:       math.MaxInt,
		Chrom:          "chr1",
		OnlyIntraChrom: true,
		OnlyInterChrom: false},
}

func TestBedPeFilter(t *testing.T) {
	var err error
	var s Settings
	for _, v := range BedpeFilterTests {
		s = Settings{
			InFile:         v.InFile,
			OutFile:        v.OutFile,
			MinScore:       v.MinScore,
			MaxScore:       v.MaxScore,
			MinDistance:    v.MinDistance,
			MaxDistance:    v.MaxDistance,
			MinStart:       v.MinStart,
			MaxStart:       v.MaxStart,
			Chrom:          v.Chrom,
			OnlyInterChrom: v.OnlyInterChrom,
			OnlyIntraChrom: v.OnlyIntraChrom,
		}
		bedpeFilter(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in bedFilter.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
