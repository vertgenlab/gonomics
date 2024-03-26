package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var HistTests = []struct {
	OboFile      string
	GafNamesFile string
	OutFile      string
	ExpectedFile string
}{
	{OboFile: "testdata/go.obo",
		GafNamesFile: "testdata/gafFiles.txt",
		OutFile:      "testdata/table.tsv",
		ExpectedFile: "testdata/expectedTable.tsv",
	},
}

func TestOboHistogram(t *testing.T) {
	for _, h := range HistTests {
		ontologyHistogram(h.OboFile, h.GafNamesFile, h.OutFile)
		if !fileio.AreEqual(h.ExpectedFile, h.OutFile) {
			t.Errorf("Error in output file for oboHistogram")
		} else {
			err := os.Remove(h.OutFile)
			exception.PanicOnErr(err)
		}
	}

}
