package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
)

var GtfToBedTests = []struct {
	inFile        string
	outFile       string
	expectedFile  string
	tss           bool
	chromSizeFile string
	merge         bool
}{
	{"testdata/test.gtf", "testdata/tmp.bed", "testdata/testOut.bed", false, "", false},
	{"testdata/test.gtf", "testdata/tmp.tss.bed", "testdata/expected.tss.bed", true, "testdata/chr1.chrom.sizes", false},
}

func TestGtfToBed(t *testing.T) {
	var err error
	for _, v := range GtfToBedTests {
		gtfToBed(v.inFile, v.outFile, v.tss, v.chromSizeFile, v.merge)
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
