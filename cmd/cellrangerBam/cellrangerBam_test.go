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
	byCell         bool
	samOut         bool
}{
	{"testdata/in.bam", "testdata/out.default.txt", "testdata/expected.default.txt", "", true, false},
	{"testdata/in.bam", "testdata/out.pseudobulk.txt", "testdata/expected.pseudobulk.txt", "", false, false},
	{"testdata/in.bam", "testdata/out.pseudobulkNormalized.txt", "testdata/expected.pseudobulkNormalized.txt", "testdata/inputNormTable.txt", false, false},
	{"testdata/in.bam", "testdata/out.samOut.sam", "testdata/expected.samOut.sam", "", false, true},
}

func TestCellrangerBam(t *testing.T) {
	var err error
	for _, v := range cellrangerBamTests {
		cellrangerBam(v.inFile, v.outFile, v.byCell, v.inputNormTable, v.samOut)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error: cellrangerBam files %s and %s are not equal to one another...", v.outFile, v.expectedFile)
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
