package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var WigStatsTests = []struct {
	InputNoGapBed    string
	ChromSizes       string
	InputWig         string
	OutputFile       string
	ExpectedFile     string
	MissingDataValue float64
}{
	{InputNoGapBed: "testdata/stats/test.noGap.bed",
		InputWig:         "testdata/stats/test.wig",
		ChromSizes:       "testdata/stats/test.chrom.sizes",
		OutputFile:       "testdata/stats/tmp.tsv",
		ExpectedFile:     "testdata/stats/expected.tsv",
		MissingDataValue: -10},
}

func TestWigStats(t *testing.T) {
	var err error
	var s StatsSettings
	for _, v := range WigStatsTests {
		s = StatsSettings{
			InFile:           v.InputWig,
			NoGapFile:        v.InputNoGapBed,
			ChromSizes:       v.ChromSizes,
			OutFile:          v.OutputFile,
			MissingDataValue: v.MissingDataValue,
		}
		wigStats(s)
		if !fileio.AreEqual(v.OutputFile, v.ExpectedFile) {
			t.Errorf("Error in wigStats.")
		} else {
			err = os.Remove(v.OutputFile)
			exception.PanicOnErr(err)
		}
	}
}
