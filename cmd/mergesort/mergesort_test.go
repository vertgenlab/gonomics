package main

import (
	"fmt"
	"os"
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var MergeSortTests = []struct {
	InFile           string
	InFileR2         string
	OutFile          string
	OutFileR2        string
	ExpectedFile     string
	ExpectedFileR2   string
	NumLinesPerChunk int
	SortCriteria     string
	fastqPE          bool
	tmpDir           string
}{
	{"testdata/bedFileTest.bed", "", "testdata/out.bed", "", "testdata/expectedSortByCoord.bed", "", 1000000, "byGenomicCoordinates", false, ""},
	{"testdata/small.sam", "", "testdata/out.sam", "", "testdata/expected.small.sam", "", 1000000, "byGenomicCoordinates", false, ""},
	{"testdata/singleCell.sam", "", "testdata/out.singleCell.sam", "", "testdata/expected.singleCell.sam", "tmp", 1000000, "singleCellBx", false, ""},
	{"testdata/test.vcf", "", "testdata/out.vcf", "", "testdata/expected.vcf", "", 1000000, "byGenomicCoordinates", false, ""},
	{"testdata/test.axt", "", "testdata/out.axt", "", "testdata/expected.axt", "", 1000000, "byGenomicCoordinates", false, ""},
	{"testdata/test_R1.fastq", "", "testdata/out_R1.fastq", "", "testdata/expected_R1.fastq", "", 1000000, "byGenomicCoordinates", false, ""},
	{"testdata/test_R1.fastq", "testdata/test_R2.fastq", "testdata/out_R1.fastq", "testdata/out_R2.fastq", "testdata/expected_R1.fastq", "testdata/expected_R2.fastq", 1000000, "byGenomicCoordinates", true, ""},
	{"testdata/small.sam", "", "testdata/out.readName.sam", "", "testdata/expected.readName.sam", "", 1000000, "readName", false, "testdata"},

	// TODO enable giraf sorting after pointers are removed
	//{"testdata/test.giraf", "testdata/out.giraf", "testdata/expected.giraf", "tmp", 1000000, "byGenomicCoordinates"},
}

func TestMergeSort(t *testing.T) {
	var err error
	for _, v := range MergeSortTests {
		if v.fastqPE {
			v.InFile = fmt.Sprintf("%s,%s", v.InFile, v.InFileR2)
			v.OutFile = fmt.Sprintf("%s,%s", v.OutFile, v.OutFileR2)
		}

		mergeSort(v.InFile, v.OutFile, v.NumLinesPerChunk, v.SortCriteria, v.tmpDir)

		if v.fastqPE {
			tmpIn := strings.Split(v.InFile, ",")
			v.InFile = tmpIn[0]
			tmpOut := strings.Split(v.OutFile, ",")
			v.OutFile = tmpOut[0]
		}
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in mergeSort: %s.", v.OutFile)
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
		if v.fastqPE {
			if !fileio.AreEqual(v.OutFileR2, v.ExpectedFileR2) {
				t.Errorf("Error in mergeSort: %s.", v.OutFileR2)
			} else {
				err = os.Remove(v.OutFileR2)
				exception.PanicOnErr(err)
			}
		}
	}
}
