package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var SamFilterTests = []struct {
	inputFile        string
	outputFile       string
	expectedFile     string
	AlignQualFilter  int
	AlignLenFilter   int
	FilterCigar      string
	CollapseByUmi    bool
	SingleCellFormat bool
	FilterByRegion   string
	FilterByFlag     int
	SortByPosition   bool
	OutBam           bool
	NoHeader         bool
}{
	{"testdata/in.test.sam", "testdata/out.noSetting.sam", "testdata/expected.noSetting.sam", 0, 0, "", false, false, "", 0, false, false, false},
	{"testdata/in.test.sam", "testdata/out.QFilter.sam", "testdata/expected.QFilter.sam", 2, 0, "", false, false, "", 0, false, false, false},
	{"testdata/in.test.sam", "testdata/out.LenFilter.sam", "testdata/expected.LenFilter.sam", 0, 70, "", false, false, "", 0, false, false, false},
	{"testdata/in.test.sam", "testdata/out.cigar.sam", "testdata/expected.cigar.sam", 0, 0, "starrSeqIntrons", false, false, "", 0, false, false, false},
	//{"testdata/in.test.sam", "testdata/out.collapseUMI.sam", "testdata/expected.collapseUMI.sam", 0, 0, "", true, true, "", 0, false, false, false},  not ready for bam function yet
	{"testdata/in.test.sam", "testdata/out.region.sam", "testdata/expected.region.sam", 0, 0, "", false, false, "chrSS", 0, false, false, false},
	{"testdata/in.test.sam", "testdata/out.flag.sam", "testdata/expected.flag.sam", 0, 0, "", false, false, "", 16, false, false, false},
	{"testdata/in.test.sam", "testdata/out.sort.sam", "testdata/expected.sort.sam", 0, 0, "", false, false, "", 0, true, false, false},
	//{"testdata/in.test.sam", "testdata/out.test.bam", "testdata/expected.test.bam", 0, 0, "", false, false, "", 0, false, true, false}, not ready for bam function yet
	{"testdata/in.test.sam", "testdata/out.noHead.sam", "testdata/expected.noHead.sam", 0, 0, "", false, false, "", 0, false, false, true},
}

func TestSamFilter(t *testing.T) {
	var err error
	var s Settings
	for _, v := range SamFilterTests {
		s = Settings{
			InFile:           v.inputFile,
			OutFile:          v.outputFile,
			AlignQualFilter:  v.AlignQualFilter,
			AlignLenFilter:   v.AlignLenFilter,
			FilterCigar:      v.FilterCigar,
			CollapseByUmi:    v.CollapseByUmi,
			SingleCellFormat: v.SingleCellFormat,
			FilterByRegion:   v.FilterByRegion,
			FilterByFlag:     v.FilterByFlag,
			SortByPosition:   v.SortByPosition,
			OutBam:           v.OutBam,
			NoHeader:         v.NoHeader,
		}
		runFilter(s)
		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error in samFilter.")
		}
		err = os.Remove(v.outputFile)
		if err != nil {
			common.ExitIfError(err)
		}
	}
}
