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
	expectedFile  string
}{
	{"testdata/testConversionTable.txt", "testdata/in.txt", "testdata/tmp.out.txt", false, false, "testdata/expected.txt"},
	{"", "testdata/in.txt", "testdata/tmp.out.txt", true, false, "testdata/expected.txt"},
}

func TestBedPeOverlap(t *testing.T) {
	var err error
	for _, v := range geneIdToNameTests {
		geneIdToName(v.geneNameTable, v.in, v.out, v.ncbi, v.ensembl)
		if !fileio.AreEqual(v.out, v.expectedFile) {
			t.Errorf("Error in geneIdToName.")
		} else {
			err = os.Remove(v.out)
			exception.PanicOnErr(err)
		}
	}
}
