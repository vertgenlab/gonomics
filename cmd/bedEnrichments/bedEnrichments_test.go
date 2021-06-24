package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedEnrichmentsTests = []struct {
	method          string
	elements1File   string
	elements2File   string
	noGapFile       string
	expectedFile    string
	trimToRefGenome bool
}{
	{"exact", "testdata/elements1.bed", "testdata/elements2.bed", "testdata/tinyNoGap.bed", "testdata/elements1.elements2.enrichment.txt", false},
	{"exact", "testdata/elements1.bed", "testdata/elements1.bed", "testdata/tinyNoGap.bed", "testdata/elements1.elements1.enrichment.txt", false},
	{"exact", "testdata/elements1.bed", "testdata/elements2.bed", "testdata/tinyNoGap.bed", "testdata/elements1.elements2.enrichment.txt", true},//should have no effect to trim to ref Genome if all elements are in the genome.
	{"exact", "testdata/elements1.bed", "testdata/elements3.bed", "testdata/tinyNoGap.bed", "testdata/elements1.elements3.enrichment.txt", true},//elements3 is elements2 with extra elements outside the genome, should be the same answer as the previous check.
}

func TestBedEnrichments(t *testing.T) {
	var err error
	for _, v := range BedEnrichmentsTests {
		bedEnrichments(v.method, v.elements1File, v.elements2File, v.noGapFile, "testdata/tmp.txt", 0, v.trimToRefGenome)
		if !fileio.AreEqual("testdata/tmp.txt", v.expectedFile) {
			t.Errorf("Error in bedEnrichments.")
		}
		err = os.Remove("testdata/tmp.txt")
		exception.PanicOnErr(err)
	}
}
