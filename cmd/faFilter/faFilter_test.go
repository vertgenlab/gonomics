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
	maxGC        float64
	minGC        float64
}{
	{"testdata/minSizeTest.fa", "testdata/minSizeOutput.fa", "testdata/minSizeExpected.fa", "", "", "", false, 0, -1, 10, 100, 0},
	{"testdata/nameContainsTest.fa", "testdata/nameContainsOutput.fa", "testdata/nameContainsExpected.fa", "", "", "_maternal", false, 0, -1, 0, 100, 0},
	{"testdata/maxGCTest.fa", "testdata/maxGCOutput.fa", "testdata/maxGCExpected.fa", "", "", "", false, 0, -1, 0, 65, 0},
	{"testdata/minGCTest.fa", "testdata/minGCOutput.fa", "testdata/minGCExpected.fa", "", "", "", false, 0, -1, 0, 100, 30},
}

func TestFaFilter(t *testing.T) {
	var err error
	for _, v := range FaFilterTests {
		faFilter(v.inputFile, v.outputFile, v.name, v.notName, v.nameContains, v.refPositions, v.start, v.end, v.minSize, v.maxGC, v.minGC)
		records := fasta.Read(v.outputFile)
		expected := fasta.Read(v.expectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error in faFilter.")
		}
		err = os.Remove(v.outputFile)
		exception.PanicOnErr(err)
	}
}
