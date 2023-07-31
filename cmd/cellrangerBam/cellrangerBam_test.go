package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var cellrangerBamTests = []struct {
	inFile         string
	outFile        string
	expectedFile   string
	inputNormTable string
	psuedobulk     bool
}{
	{"testdata/in.bam", "testdata/out.default.txt", "testdata/expected.default.txt", "", false},
	{"testdata/in.bam", "testdata/out.pseudobulk.txt", "testdata/expected.pseudobulk.txt", "", true},
	{"testdata/in.bam", "testdata/out.pseudobulkNormalized.txt", "testdata/expected.pseudobulkNormalized.txt", "testdata/inputNormTable.txt", true},
}

func TestCellrangerBam(t *testing.T) {
	var err error
	for _, v := range cellrangerBamTests {
		cellrangerBam(v.inFile, v.outFile, v.psuedobulk, v.inputNormTable)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error: cellrangerBam files %s and %s are not equal to one another...", v.outFile, v.expectedFile)
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
