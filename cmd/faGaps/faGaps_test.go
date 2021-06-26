package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedOutTests = []struct {
	inFile       string
	expectedFile string
}{
	{"testdata/testBedOut.fa", "testdata/testBedOut.noGap.bed"},
}

func TestFaGapsBedOut(t *testing.T) {
	var err error
	for _, v := range BedOutTests {
		faNoGap(v.inFile, "testdata/tmp.txt")
		if !fileio.AreEqual(v.expectedFile, "testdata/tmp.txt") {
			t.Errorf("Error in FaGaps NoGapBed option.")
		} else {
			err = os.Remove("testdata/tmp.txt")
			exception.PanicOnErr(err)
		}
	}
}
