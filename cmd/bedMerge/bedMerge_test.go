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
	MergeAdjacent bool
	LowMem bool
}{
	{"testdata/test.bed", "testdata/test.merged.bed", false, false},
	{"testdata/test.bed", "testdata/test.adjacent.merged.bed", true, false},
	{"testdata/test.presorted.bed", "testdata/test.lowmem.merged.bed", false, true},
	{"testdata/test.presorted.bed", "testdata/test.adjacent.lowmem.merged.bed", true, true},
}

func TestBedMerge(t *testing.T) {
	var err error
	for _, i := range BedMergeTests {
		bedMerge(i.InFile, "testdata/tmp.txt", i.MergeAdjacent, i.LowMem)
		if !fileio.AreEqual("testdata/tmp.txt", i.ExpectedFile) {
			t.Errorf("Error in bedMerge.")
		} else {
			err = os.Remove("testdata/tmp.txt")
			exception.PanicOnErr(err)
		}
	}
}
