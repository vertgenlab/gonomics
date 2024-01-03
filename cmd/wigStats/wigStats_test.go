package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var WigStatsTests = []struct {
	InputNoGapBed    string
	ChromSizes       string
	InputWig         string
	OutputFile       string
	ExpectedFile     string
	MissingDataValue float64
}{
	{InputNoGapBed: "testdata/test.noGap.bed",
		InputWig:         "testdata/test.wig",
		ChromSizes:       "testdata/test.chrom.sizes",
		OutputFile:       "testdata/tmp.tsv",
		ExpectedFile:     "testdata/expected.tsv",
		MissingDataValue: -10},
}

func TestWigStats(t *testing.T) {
	var err error
	var s Settings
	for _, v := range WigStatsTests {
		s = Settings{
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
