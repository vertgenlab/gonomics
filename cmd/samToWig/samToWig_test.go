package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	//"github.com/vertgenlab/gonomics/bed" //will only use when bed.AllAreEqual can accomodate comparing not just the first 3 fields in the future
	"os"
	"testing"
)

var SamToWigTests = []struct {
	inFile           string
	reference	string
	outFile_expected string
	fragLength       int
}{
	{"testdata/test1.sam", "testdata/test.chrom.sizes", "testdata/test1.wig", -1},
	{"testdata/test2.sam", "testdata/test.chrom.sizes", "testdata/test2.wig", 30},
}

func TestSamToWig(t *testing.T) {
	for _, v := range SamToWigTests {
		samToWig(v.inFile, v.reference, "outFile_tmp.wig", v.fragLength)
		if !fileio.AreEqual("outFile_tmp.wig", v.outFile_expected) {
			t.Errorf("Error in samToWig")
		}
		err := os.Remove("outFile_tmp.wig")
		if err != nil {
			common.ExitIfError(err)
		}
	}
}
