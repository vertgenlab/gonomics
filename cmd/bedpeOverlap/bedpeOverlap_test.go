package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var BedPeOverlapTests = []struct {
	selectFile       string
	inBedPe          string
	outBedPe         string
	bedSelect        bool
	overlapBoth      bool
	overlapThreshold float64
	keepNames        bool
	expectedFile     string
}{
	{"testdata/selectBedPe.bedpe", "testdata/inBedPe.bedpe", "testdata/tmp.bedpe", false, false, 0, false, "testdata/expected.bedpe"},
	{"testdata/select.bed", "testdata/inBedPe.bedpe", "testdata/tmp.bedSelect.bedpe", true, false, 0, false, "testdata/expected.bedSelect.bedpe"},
	{"testdata/select.names.bed", "testdata/inBedPe.bedpe", "testdata/tmp.bedSelect.names.bedpe", true, false, 0, true, "testdata/expected.bedSelect.names.bedpe"},
	{"testdata/select.bed", "testdata/inBedPe.bedpe", "testdata/tmp.bedSelectOverlapThresh.bedpe", true, false, 0.5, false, "testdata/expected.bedSelect.overlapThresh.bedpe"},
	{"testdata/selectBedBoth.bed", "testdata/inBedPe.bedpe", "testdata/tmp.bedSelectOverlapBoth.bedpe", true, true, 0, false, "testdata/expected.bedSelect.both.bedpe"},
	{"testdata/selectBedBothThresh.bed", "testdata/inBedPe.bedpe", "testdata/tmp.bedSelectOverlapThreshOverlapBoth.bedpe", true, true, 0.5, false, "testdata/expected.bedSelect.both.bedpe"},
}

func TestBedPeOverlap(t *testing.T) {
	var err error
	for _, v := range BedPeOverlapTests {
		bedpeOverlap(v.selectFile, v.inBedPe, v.outBedPe, v.bedSelect, v.overlapThreshold, v.overlapBoth, v.keepNames)
		if !fileio.AreEqual(v.outBedPe, v.expectedFile) {
			t.Errorf("Error: bedpeOverlap files %s and %s are not equal to one another...", v.outBedPe, v.expectedFile)
		} else {
			err = os.Remove(v.outBedPe)
			exception.PanicOnErr(err)
		}
	}
}
