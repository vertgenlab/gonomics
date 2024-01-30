package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/fasta"
)

var FaFilterTests = []struct {
	inputFile    string
	outputFile   string
	expectedFile string
	name         string
	notName      string
	nameContains string
	refPositions bool
	start        int
	end          int
	minSize      int
	maxGC        int
}{
	{"testdata/minSizeTest.fa", "testdata/minSizeOutput.fa", "testdata/minSizeExpected.fa", "", "", "", false, 0, -1, 10, 100},
	{"testdata/nameContainsTest.fa", "testdata/nameContainsOutput.fa", "testdata/nameContainsExpected.fa", "", "", "_maternal", false, 0, -1, 0, 100},
	{"testdata/maxGCTest.fa", "testdata/maxGCOutput.fa", "testdata/maxGCExpected.fa", "", "", "", false, 0, -1, 0, 65},
}

func TestFaFilter(t *testing.T) {
	var err error
	for _, v := range FaFilterTests {
		faFilter(v.inputFile, v.outputFile, v.name, v.notName, v.nameContains, v.refPositions, v.start, v.end, v.minSize, v.maxGC)
		records := fasta.Read(v.outputFile)
		expected := fasta.Read(v.expectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error in faFilter.")
		}
		err = os.Remove(v.outputFile)
		exception.PanicOnErr(err)
	}
}
