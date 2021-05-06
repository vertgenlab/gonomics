package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"os"
	"testing"
)

var MultiFaExtractTests = []struct {
	inputFile string
	outputFile string
	expectedFile string
	removeGaps bool
	start int
	end int
}{
	{"testdata/testInput.fa", "testdata/testOut.fa", "testdata/testOut.10to200.fa", false, 10, 200},
	{"testdata/testInput.fa", "testdata/testOut.fa", "testdata/testOut.10to200.RemoveGaps.fa", true, 10, 200},
}

func TestMultiFaExtract(t *testing.T) {
	var err error
	for _, v := range MultiFaExtractTests {
		multiFaExtract(v.inputFile, v.outputFile, v.start, v.end, v.removeGaps)
		records := fasta.Read(v.outputFile)
		expected := fasta.Read(v.expectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error in multiFaExtract.")
		}
		err = os.Remove(v.outputFile)
		if err != nil {
			common.ExitIfError(err)
		}
	}
}