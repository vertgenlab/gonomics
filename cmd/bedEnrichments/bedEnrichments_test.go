package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedEnrichmentsTests = []struct {
	method string
	elements1File     string
	elements2File     string
	noGapFile         string
	expectedFile      string
}{
	{"exact","testdata/elements1.bed", "testdata/elements2.bed", "testdata/tinyNoGap.bed", "testdata/elements1.elements2.enrichment.txt"},
	{"exact","testdata/elements1.bed", "testdata/elements1.bed", "testdata/tinyNoGap.bed", "testdata/elements1.elements1.enrichment.txt"},
}

func TestBedEnrichments(t *testing.T) {
	var err error
	for _, v := range BedEnrichmentsTests {
		bedEnrichments(v.method, v.elements1File, v.elements2File, v.noGapFile, "testdata/tmp.txt", 0)
		if !fileio.AreEqual("testdata/tmp.txt", v.expectedFile) {
			t.Errorf("Error in bedEnrichments.")
		}
		err = os.Remove("testdata/tmp.txt")
		exception.PanicOnErr(err)
	}
}
