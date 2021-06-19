package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedEnrichmentsTests = []struct {
	elements1File string
	elements2File string
	noGapFile     string
	expectedFile  string
	approximateMethod bool
}{
	{"testdata/elements1.bed", "testdata/elements2.bed", "testdata/tinyNoGap.bed", "testdata/elements1.elements2.enrichment.txt", false},
	{"testdata/elements1.bed", "testdata/elements1.bed", "testdata/tinyNoGap.bed", "testdata/elements1.elements1.enrichment.txt", false},
}

func TestBedEnrichments(t *testing.T) {
	var err error
	for _, v := range BedEnrichmentsTests {
		bedEnrichments(v.elements1File, v.elements2File, v.noGapFile, "testdata/tmp.txt", v.approximateMethod)
		if !fileio.AreEqual("testdata/tmp.txt", v.expectedFile) {
			t.Errorf("Error in bedEnrichments.")
		}
		err = os.Remove("testdata/tmp.txt")
		exception.PanicOnErr(err)
	}
}
