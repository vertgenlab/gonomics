package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var MergeSortTests = []struct {
	InFile           string
	OutFile          string
	ExpectedFile     string
	TmpFilePrefix    string
	NumLinesPerChunk int
	SortCriteria     string
}{
	{"testdata/bedFileTest.bed", "testdata/out.bed", "testdata/expectedSortByCoord.bed", "tmp", 1000000, "byGenomicCoordinates"},
	{"testdata/small.sam", "testdata/out.sam", "testdata/expected.small.sam", "tmp", 1000000, "byGenomicCoordinates"},
	{"testdata/singleCell.sam", "testdata/out.singleCell.sam", "testdata/expected.singleCell.sam", "tmp", 1000000, "singleCellBx"},
	{"testdata/test.vcf", "testdata/out.vcf", "testdata/expected.vcf", "tmp", 1000000, "byGenomicCoordinates"},
	{"testdata/test.axt", "testdata/out.axt", "testdata/expected.axt", "tmp", 1000000, "byGenomicCoordinates"},
	// TODO enable giraf sorting after pointers are removed
	//{"testdata/test.giraf", "testdata/out.giraf", "testdata/expected.giraf", "tmp", 1000000, "byGenomicCoordinates"},
}

func TestMergeSort(t *testing.T) {
	var err error
	for _, v := range MergeSortTests {
		mergeSort(v.InFile, v.OutFile, v.NumLinesPerChunk, v.SortCriteria)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in mergeSort: %s.", v.OutFile)
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
