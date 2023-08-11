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
	singleCell     string
	binCells       int
}{
	{"testdata/in.bam", "testdata/out.default.txt", "testdata/expected.default.txt", "", true, false, "", 0},
	{"testdata/in.bam", "testdata/out.pseudobulk.txt", "testdata/expected.pseudobulk.txt", "", false, false, "", 0},
	{"testdata/in.bam", "testdata/out.pseudobulkNormalized.txt", "testdata/expected.pseudobulkNormalized.txt", "testdata/inputNormTable.txt", false, false, "", 0},
	//{"testdata/in.bam", "testdata/out.samOut.sam", "testdata/expected.samOut.sam", "", false, true, "", 0},
	{"testdata/in.bam", "testdata/out.singleCell.txt", "testdata/expected.singleCell.txt", "", false, false, "testdata/cellTypes.txt", 0},
	{"testdata/in.bam", "testdata/out.singleCellNormalized.txt", "testdata/expected.singleCellNormalized.txt", "testdata/inputNormTable.txt", false, false, "testdata/cellTypes.txt", 0},
}

func TestCellrangerBam(t *testing.T) {
	var err error
	for _, v := range cellrangerBamTests {
		parseBam(v.inFile, v.outFile, v.byCell, v.inputNormTable, v.samOut, v.singleCell, v.binCells)
		if !fileio.AreEqual(v.outFile, v.expectedFile) {
			t.Errorf("Error: cellrangerBam files %s and %s are not equal to one another...", v.outFile, v.expectedFile)
		} else {
			err = os.Remove(v.outFile)
			exception.PanicOnErr(err)
		}
	}
}
