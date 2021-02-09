package main

import (
	"testing"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"os"
)

var GtfToBedTests = []struct {
	inFile	string
	expectedFile	string
}{
	{"testdata/test.gtf", "testdata/testOut.bed"},
}

func TestGtfToBed(t *testing.T) {
	for _, v := range GtfToBedTests {
		gtfToBed(v.inFile, "tmp.bed")
		records := bed.Read("tmp.bed")
		expected := bed.Read(v.expectedFile)
		if !bed.AllAreEqual(records, expected) {
			t.Errorf("Error in gtfToBed.")
		}
	}
	err := os.Remove("tmp.bed")
	if err != nil {
		common.ExitIfError(err)
	}
}
