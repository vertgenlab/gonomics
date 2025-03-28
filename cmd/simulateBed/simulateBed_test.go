package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var SimulateBedTests = []struct {
	RegionCount  int
	SimLength    int
	MatchedBed   string
	NoGapFile    string
	OutFile      string
	SetSeed      int64
	ExpectedFile string
}{
	{10, 1000, "", "testdata/test.noGap.bed", "testdata/tmp.bed", 10, "testdata/expected.bed"},
	{0, 0, "testdata/expected.bed", "testdata/test.noGap.bed", "testdata/tmp.bed", 10, "testdata/expected.matched.bed"},
}

func TestSimulateBed(t *testing.T) {
	var err error
	for _, v := range SimulateBedTests {
		simulateBed(v.RegionCount, v.SimLength, v.MatchedBed, v.NoGapFile, v.OutFile, v.SetSeed)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in SimulateBed. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
