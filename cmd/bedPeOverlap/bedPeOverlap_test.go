package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedPeOverlapTests = []struct {
	selectFile   string
	inBedPe      string
	outBedPe     string
	bedSelect    bool
	expectedFile string
}{
	{"testdata/selectBedPe.bedpe", "testdata/inBedPe.bedpe", "testdata/tmp.bedpe", false, "testdata/expected.bedpe"},
	{"testdata/select.bed", "testdata/inBedPe.bedpe", "testdata/tmp.bedSelect.bedpe", true, "testdata/expected.bedSelect.bedpe"},
}

func TestBedPeOverlap(t *testing.T) {
	var err error
	for _, v := range BedPeOverlapTests {
		bedpeOverlap(v.selectFile, v.inBedPe, v.outBedPe, v.bedSelect)
		if !fileio.AreEqual(v.outBedPe, v.expectedFile) {
			t.Errorf("Error in bedPeOverlap.")
		} else {
			err = os.Remove(v.outBedPe)
			exception.PanicOnErr(err)
		}
	}
}
