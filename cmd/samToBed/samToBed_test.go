package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var SamToBedTests = []struct {
	inFile           string
	outFile_expected string
	fragLength       int
}{
	{"testdata/test1.sam", "testdata/test1.bed", -1},
	{"testdata/test2.sam", "testdata/test2.bed", 30},
}

func TestSamToBed(t *testing.T) {
	for _, v := range SamToBedTests {
		samToBed(v.inFile, "outFile_tmp.bed", v.fragLength)
		if !fileio.AreEqual("outFile_tmp.bed", v.outFile_expected) {
			t.Errorf("Error in samToBed")
		}
		err := os.Remove("outFile_tmp.bed")
		exception.PanicOnErr(err)
	}
}
