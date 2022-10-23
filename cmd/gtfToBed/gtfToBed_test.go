package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"testing"
)

var GtfToBedTests = []struct {
	inFile        string
	outFile       string
	expectedFile  string
	tss           bool
	chromSizeFile string
}{
	{"testdata/test.gtf", "testdata/tmp.bed", "testdata/testOut.bed", false, ""},
	{"testdata/test.gtf", "testdata/tmp.tss.bed", "testdata/expected.tss.bed", true, "testdata/chr1.chrom.sizes"},
}

func TestGtfToBed(t *testing.T) {
	var err error
	for _, v := range GtfToBedTests {
		gtfToBed(v.inFile, v.outFile, v.tss, v.chromSizeFile)
		records := bed.Read(v.outFile)
		expected := bed.Read(v.expectedFile)
		if !bed.AllAreEqual(records, expected) {
			t.Errorf("Error in gtfToBed.")
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
