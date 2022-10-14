package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedPeOverlapTests = []struct {
	selectBedPe  string
	inBedPe      string
	outBedPe     string
	expectedFile string
}{
	{"testdata/selectBedPe.bedpe", "testdata/inBedPe.bedpe", "testdata/tmp.bedpe", "testdata/expected.bedpe"},
}

func TestBedPeOverlap(t *testing.T) {
	var err error
	for _, v := range BedPeOverlapTests {
		bedpeOverlap(v.selectBedPe, v.inBedPe, v.outBedPe)
		if !fileio.AreEqual(v.outBedPe, v.expectedFile) {
			t.Errorf("Error in bedPeOverlap.")
		} else {
			err = os.Remove(v.outBedPe)
			exception.PanicOnErr(err)
		}
	}
}
