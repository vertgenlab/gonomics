package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedMergeTests = []struct {
	InFile       string
	ExpectedFile string
	MergeAdjacent bool
}{
	{"testdata/test.bed", "testdata/test.merged.bed", false},
	{"testdata/test.bed", "testdata/test.adjacent.merged.bed", true},
}

func TestBedMerge(t *testing.T) {
	var err error
	for _, i := range BedMergeTests {
		bedMerge(i.InFile, "tmp.txt", i.MergeAdjacent)
		if !fileio.AreEqual("tmp.txt", i.ExpectedFile) {
			t.Errorf("Error in bedMerge.")
		} else {
			err = os.Remove("tmp.txt")
			exception.PanicOnErr(err)
		}
	}
}
