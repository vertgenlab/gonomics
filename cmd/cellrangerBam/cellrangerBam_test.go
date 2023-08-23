package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var cellrangerBamTests = []struct {
	InFile         string
	OutFile        string
	ExpectedFile   string
	InputNormTable string
	ByCell         bool
	SamOut         bool
	SingleCell     string
	BinCells       int
	UmiSat         bool
	GfpNorm        string
	BedRegion      string
}{
	{"testdata/in.bam", "testdata/out.default.txt", "testdata/expected.default.txt", "", true, false, "", 0, false, "", ""},
	{"testdata/in.bam", "testdata/out.pseudobulk.txt", "testdata/expected.pseudobulk.txt", "", false, false, "", 0, false, "", ""},
	{"testdata/in.bam", "testdata/out.pseudobulkNormalized.txt", "testdata/expected.pseudobulkNormalized.txt", "testdata/inputNormTable.txt", false, false, "", 0, false, "", ""},
	//{"testdata/in.bam", "testdata/out.samOut.sam", "testdata/expected.samOut.sam", "", false, true, "", 0, false},
	{"testdata/in.bam", "testdata/out.singleCell.txt", "testdata/expected.singleCell.txt", "", false, false, "testdata/cellTypes.txt", 0, false, "", ""},
	{"testdata/in.bam", "testdata/out.singleCellNormalized.txt", "testdata/expected.singleCellNormalized.txt", "testdata/inputNormTable.txt", false, false, "testdata/cellTypes.txt", 0, false, "", ""},
	{"testdata/in.bam", "testdata/out.singleCellGfpNormalized.txt", "testdata/expected.singleCellGfpNormalized.txt", "testdata/inputNormTable.txt", false, false, "testdata/cellTypes.txt", 0, false, "testdata/in.gfp.bam", ""},
	{"testdata/in.bam", "testdata/out.pseudobulkNormalized.txt", "testdata/expected.pseudobulkNormalized.txt", "testdata/inputNormTable.txt", false, false, "", 0, false, "", "testdata/in.bed"},
}

func TestCellrangerBam(t *testing.T) {
	var err error
	var s Settings
	for _, v := range cellrangerBamTests {
		s = Settings{
			inFile:     v.InFile,
			outFile:    v.OutFile,
			normalize:  v.InputNormTable,
			byCell:     v.ByCell,
			samOut:     v.SamOut,
			scAnalysis: v.SingleCell,
			binCells:   v.BinCells,
			umiSat:     v.UmiSat,
			gfpNorm:    v.GfpNorm,
			bed:        v.BedRegion,
		}
		parseBam(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: cellrangerBam files %s and %s are not equal to one another...", v.OutFile, v.ExpectedFile)
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
