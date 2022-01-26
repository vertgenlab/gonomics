package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var SimulateBedTests = []struct {
	RegionCount  int
	SimLength    int
	NoGapFile    string
	OutFile      string
	SetSeed      int64
	ExpectedFile string
}{
	{10, 1000, "testdata/test.noGap.bed", "testdata/tmp.bed", 10, "testdata/expected.bed"},
}

func TestSimulateBed(t *testing.T) {
	var err error
	for _, v := range SimulateBedTests {
		simulateBed(v.RegionCount, v.SimLength, v.NoGapFile, v.OutFile, v.SetSeed)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in SimulateBed. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
