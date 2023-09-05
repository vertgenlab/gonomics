package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var WigStatsTests = []struct {
	inputNoGapBed    string
	inputWig         string
	outputFile       string
	expectedFile     string
	missingDataValue float64
}{
	{"testdata/test.noGap.bed", "testdata/test.wig", "testdata/tmp.tsv", "testdata/expected.tsv", -10},
}

func TestWigStats(t *testing.T) {
	var err error
	for _, v := range WigStatsTests {
		wigStats(v.inputWig, v.inputNoGapBed, v.outputFile, v.missingDataValue)
		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error in wigStats.")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
	}
}
