package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedMergeTests = []struct {
	InFile        string
	ExpectedFile  string
	OutFile       string
	MergeAdjacent bool
	Pad           int
	LowMem        bool
	KeepAllNames  bool
}{
	{"testdata/test.bed", "testdata/test.merged.bed", "testdata/out.merged.bed", false, -1, false, false},
	{"testdata/test.bed", "testdata/test.adjacent.merged.bed", "testdata/out.adjacent.mereged.bed", true, -1, false, false},
	{"testdata/test.presorted.bed", "testdata/test.lowmem.merged.bed", "testdata/out.lowmem.merged.bed", false, -1, true, false},
	{"testdata/test.presorted.bed", "testdata/test.adjacent.lowmem.merged.bed", "testdata/out.adjacent.lowmem.merged.bed", true, -1, true, false},
	{"testdata/test.names.bed", "testdata/test.names.merged.bed", "testdata/out.names.merged.bed", false, -1, false, true},
	{"testdata/test.names.bed", "testdata/test.names.adjacent.merged.bed", "testdata/out.names.adjacent.merged.bed", true, -1, false, true},
	{"testdata/testPad.presorted.bed", "testdata/test.pad.merged.bed", "testdata/out.pad.merged.bed", true, 5, true, false},
	{"testdata/testPad.presorted.bed", "testdata/test.names.pad.merged.bed", "testdata/out.names.pad.merged.bed", true, 5, false, true},
}

func TestBedMerge(t *testing.T) {
	var err error
	var dist int = -1
	for _, i := range BedMergeTests {
		dist = -1
		if i.Pad > -1 {
			dist = i.Pad + 1
		}
		if i.MergeAdjacent && i.Pad < 0 {
			dist = 1
		}
		bedMerge(i.InFile, i.OutFile, dist, i.LowMem, i.KeepAllNames)
		if !fileio.AreEqual(i.OutFile, i.ExpectedFile) {
			t.Errorf("Error in bedMerge.")
		} else {
			err = os.Remove(i.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
