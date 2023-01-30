package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedPeOverlapTests = []struct {
	selectFile       string
	inBedPe          string
	outBedPe         string
	bedSelect        bool
	overlapThreshold float64
	expectedFile     string
}{
	{"testdata/selectBedPe.bedpe", "testdata/inBedPe.bedpe", "testdata/tmp.bedpe", false, 0, "testdata/expected.bedpe"},
	{"testdata/select.bed", "testdata/inBedPe.bedpe", "testdata/tmp.bedSelect.bedpe", true, 0, "testdata/expected.bedSelect.bedpe"},
	{"testdata/select.bed", "testdata/inBedPe.bedpe", "testdata/tmp.bedSelectOverlapThresh.bedpe", true, 0.5, "testdata/expected.bedSelect.overlapThresh.bedpe"},
}

func TestBedPeOverlap(t *testing.T) {
	var err error
	for _, v := range BedPeOverlapTests {
		bedpeOverlap(v.selectFile, v.inBedPe, v.outBedPe, v.bedSelect, v.overlapThreshold)
		if !fileio.AreEqual(v.outBedPe, v.expectedFile) {
			t.Errorf("Error in bedPeOverlap.")
		} else {
			err = os.Remove(v.outBedPe)
			exception.PanicOnErr(err)
		}
	}
}
