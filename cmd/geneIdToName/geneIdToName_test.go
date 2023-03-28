package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var geneIdToNameTests = []struct {
	geneNameTable string
	in            string
	out           string
	ncbi          bool
	ensembl       bool
	keepMatching  bool
	expectedFile  string
}{
	{"testdata/testConversionTable.txt", "testdata/in.txt", "testdata/tmp.out.txt", false, false, false, "testdata/expected.txt"},
	{"", "testdata/in.txt", "testdata/tmp.out.txt", true, false, false, "testdata/expected.txt"},
	{"", "testdata/in.txt", "testdata/tmp.keepMatching.out.txt", true, false, true, "testdata/expectedKeepMatching.txt"},
}

func TestBedPeOverlap(t *testing.T) {
	var err error
	for _, v := range geneIdToNameTests {
		geneIdToName(v.geneNameTable, v.in, v.out, v.ncbi, v.ensembl, v.keepMatching)
		if !fileio.AreEqual(v.out, v.expectedFile) {
			t.Errorf("Error in geneIdToName.")
		} else {
			err = os.Remove(v.out)
			exception.PanicOnErr(err)
		}
	}
}
