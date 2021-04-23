package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"os"
	"testing"
)

var FaFilterTests = []struct {
	inputFile				string
	outputFile				string
	expectedFile			string
	name					string
	notName					string
	refPositions			bool
	start 					int
	end 					int
	minSize					int
}{
	{"testdata/minSizeTest.fa","testdata/minSizeOutput.fa","testdata/minSizeExpected.fa", "","",false,0,-1,10 },
}

func TestFaFilter(t *testing.T) {
	var err error
	for _, v := range FaFilterTests {
		faFilter(v.inputFile, v.outputFile, v.name, v.notName, v.refPositions, v.start, v.end, v.minSize)
		records := fasta.Read(v.outputFile)
		expected := fasta.Read(v.expectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error in faFilter.")
		}
		err = os.Remove(v.outputFile)
		if err != nil {
			common.ExitIfError(err)
		}
	}
}