package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/starrSeq"
	"os"
	"testing"
)

var cellrangerBamTests = []struct {
	InFile         string
	OutFile        string
	ExpectedFile   string
	InputNormTable string
	ByCell         string
	SamOut         string
	SingleCell     string
	BinCells       int
	UmiSat         string
	TransNorm      string
	BedRegion      string
	NoOut          bool
	ScCount        string
	AltMapping     string
}{
	{"testdata/in.bam", "", "testdata/expected.default.txt", "", "testdata/out.default.txt", "", "", 0, "", "", "", true, "", ""},
	{"testdata/in.bam", "testdata/out.pseudobulk.txt", "testdata/expected.pseudobulk.txt", "", "", "", "", 0, "", "", "", false, "", ""},
	{"testdata/in.bam", "testdata/out.pseudobulkNormalized.txt", "testdata/expected.pseudobulkNormalized.txt", "testdata/inputNormTable.txt", "", "false", "", 0, "", "", "", false, "", ""},
	//{"testdata/in.bam", "testdata/out.samOut.sam", "testdata/expected.samOut.sam", "", false, true, "", 0, false},
	{"testdata/in.bam", "testdata/out.singleCell.txt", "testdata/expected.singleCell.txt", "", "", "", "testdata/cellTypes.txt", 0, "", "", "", false, "", ""},
	{"testdata/in.bam", "testdata/out.singleCellNormalized.txt", "testdata/expected.singleCellNormalized.txt", "testdata/inputNormTable.txt", "", "", "testdata/cellTypes.txt", 0, "", "", "", false, "", ""},
	{"testdata/in.bam", "testdata/out.singleCellGfpNormalized.txt", "testdata/expected.singleCellGfpNormalized.txt", "testdata/inputNormTable.txt", "", "", "testdata/cellTypes.txt", 0, "", "testdata/in.gfp.bam", "", false, "", ""},
	{"testdata/in.bam", "testdata/out.pseudobulkNormalized.txt", "testdata/expected.pseudobulkNormalized.txt", "testdata/inputNormTable.txt", "", "", "", 0, "", "", "testdata/in.bed", false, "", ""},
	{"testdata/in.bam", "", "testdata/expected.countMatrix.txt", "", "", "", "", 0, "", "", "", true, "testdata/out.countMatrix.txt", ""},
	{"testdata/in.bam", "", "testdata/expected.countMatrixNorm.txt", "testdata/inputNormTable.txt", "", "", "", 0, "", "", "", true, "testdata/out.countMatrixNorm.txt", ""},
	{"testdata/in.bam", "", "testdata/expected.countMatrixGFP_cellType.txt", "", "", "", "testdata/cellTypes.txt", 0, "", "testdata/in.gfp.bam", "", true, "testdata/out.countMatrixNorm.txt", ""},
	{"testdata/in.bam", "testdata/out.altMapping.txt", "testdata/expected.altMapping.txt", "testdata/inputNormTable.txt", "", "", "", 0, "", "", "", false, "", "testdata/in.bed"},
}

func TestCellrangerBam(t *testing.T) {
	var err error
	var s starrSeq.ScStarrSeqSettings
	for _, v := range cellrangerBamTests {
		s = starrSeq.ScStarrSeqSettings{
			InFile:           v.InFile,
			OutFile:          v.OutFile,
			InputNormalize:   v.InputNormTable,
			ByCell:           v.ByCell,
			SamOut:           v.SamOut,
			ScAnalysis:       v.SingleCell,
			BinCells:         v.BinCells,
			UmiSat:           v.UmiSat,
			TransfectionNorm: v.TransNorm,
			Bed:              v.BedRegion,
			NoOut:            v.NoOut,
			CountMatrix:      v.ScCount,
			AltMapping:       v.AltMapping,
		}
		if v.AltMapping == "" {
			parseBam(s)
		} else {
			starrSeq.Alt(s)
		}
		if v.OutFile != "" {
			if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
				t.Errorf("Error: cellrangerBam files %s and %s are not equal to one another...", v.OutFile, v.ExpectedFile)
			} else {
				err = os.Remove(v.OutFile)
				exception.PanicOnErr(err)
			}
		} else {
			switch {
			case v.SingleCell != "" && v.ScCount != "" && v.TransNorm != "":
				if !fileio.AreEqual(v.ScCount, v.ExpectedFile) {
					t.Errorf("Error: cellrangerBam files %s and %s are not equal to one another...", v.ScCount, v.ExpectedFile)
				} else {
					err = os.Remove(v.ScCount)
					exception.PanicOnErr(err)
				}
			case v.ScCount != "" && v.InputNormTable == "":
				if !fileio.AreEqual(v.ScCount, v.ExpectedFile) {
					t.Errorf("Error: cellrangerBam files %s and %s are not equal to one another...", v.ScCount, v.ExpectedFile)
				} else {
					err = os.Remove(v.ScCount)
					exception.PanicOnErr(err)
				}
			case v.ScCount != "" && v.InputNormTable != "":
				if !fileio.AreEqual(v.ScCount, v.ExpectedFile) {
					t.Errorf("Error: cellrangerBam files %s and %s are not equal to one another...", v.ScCount, v.ExpectedFile)
				} else {
					err = os.Remove(v.ScCount)
					exception.PanicOnErr(err)
				}
			case v.ByCell != "":
				if !fileio.AreEqual(v.ByCell, v.ExpectedFile) {
					t.Errorf("Error: cellrangerBam files %s and %s are not equal to one another...", v.ByCell, v.ExpectedFile)
				} else {
					err = os.Remove(v.ByCell)
					exception.PanicOnErr(err)
				}
			}
		}
	}
}
