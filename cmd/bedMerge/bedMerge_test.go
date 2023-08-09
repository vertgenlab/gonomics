package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var BedMergeTests = []struct {
	InFile        string
	ExpectedFile  string
	MergeAdjacent bool
	LowMem        bool
	KeepAllNames  bool
}{
	{"testdata/test.bed", "testdata/test.merged.bed", false, false, false},
	{"testdata/test.bed", "testdata/test.adjacent.merged.bed", true, false, false},
	{"testdata/test.presorted.bed", "testdata/test.lowmem.merged.bed", false, true, false},
	{"testdata/test.presorted.bed", "testdata/test.adjacent.lowmem.merged.bed", true, true, false},
	{"testdata/test.names.bed", "testdata/test.names.merged.bed", false, false, true},
	{"testdata/test.names.bed", "testdata/test.names.adjacent.merged.bed", true, false, true},
}

func TestBedMerge(t *testing.T) {
	var err error
	for _, i := range BedMergeTests {
		bedMerge(i.InFile, "testdata/tmp.txt", i.MergeAdjacent, i.LowMem, i.KeepAllNames)
		if !fileio.AreEqual("testdata/tmp.txt", i.ExpectedFile) {
			t.Errorf("Error in bedMerge.")
		} else {
			err = os.Remove("testdata/tmp.txt")
			exception.PanicOnErr(err)
		}
	}
}
