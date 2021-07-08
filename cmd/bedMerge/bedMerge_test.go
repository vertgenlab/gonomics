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
}{
	{"testdata/test.bed", "testdata/test.merged.bed"},
}

func TestBedMerge(t *testing.T) {
	var err error
	for _, i := range BedMergeTests {
		bedMerge(i.InFile, "tmp.txt")
		if !fileio.AreEqual("tmp.txt", i.ExpectedFile) {
			t.Errorf("Error in bedMerge.")
		} else {
			err = os.Remove("tmp.txt")
			exception.PanicOnErr(err)
		}
	}
}
