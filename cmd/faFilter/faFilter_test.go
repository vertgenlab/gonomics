package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/fasta"
)

var FaFilterTests = []struct {
	inputFile      string
	outputFile     string
	expectedFile   string
	name           string
	notName        string
	nameContains   string
	refPositions   bool
	start          int
	end            int
	minSize        int
	maxGC          float64
	minGC          float64
	finalBases     int
	cutFinalNbases int
	append         string
}{
	{"testdata/minSizeTest.fa", "testdata/minSizeOutput.fa", "testdata/minSizeExpected.fa", "", "", "", false, 0, -1, 10, 100, 0, -1, -1, ""},
	{"testdata/nameContainsTest.fa", "testdata/nameContainsOutput.fa", "testdata/nameContainsExpected.fa", "", "", "_maternal", false, 0, -1, 0, 100, 0, -1, -1, ""},
	{"testdata/maxGCTest.fa", "testdata/maxGCOutput.fa", "testdata/maxGCExpected.fa", "", "", "", false, 0, -1, 0, 65, 0, -1, -1, ""},
	{"testdata/minGCTest.fa", "testdata/minGCOutput.fa", "testdata/minGCExpected.fa", "", "", "", false, 0, -1, 0, 100, 30, -1, -1, ""},
	{"testdata/nameContainsTest.fa", "testdata/finalNbasesOut.fa", "testdata/finalNbasesExpected.fa", "", "", "", false, 0, -1, 0, 100, 0, 5, -1, ""},
	{"testdata/nameContainsTest.fa", "testdata/cutFinalNbasesOut.fa", "testdata/cutFinalNbasesExpected.fa", "", "", "", false, 0, -1, 0, 100, 0, -1, 5, ""},
	{"testdata/minSizeTest.fa", "testdata/appendOutput.fa", "testdata/appendExpected.fa", "", "", "", false, 0, -1, 0, 100, 0, -1, -1, "testdata/appendSeq.fa"},
}

func TestFaFilter(t *testing.T) {
	var err error
	var s settings
	for _, v := range FaFilterTests {
		s = settings{
			InFile:         v.inputFile,
			OutFile:        v.outputFile,
			Name:           v.name,
			NotName:        v.notName,
			NameContains:   v.nameContains,
			MaxGC:          v.maxGC,
			MinGC:          v.minGC,
			CutFinalNBases: v.cutFinalNbases,
			FinalNBases:    v.finalBases,
			Start:          v.start,
			End:            v.end,
			RefPositions:   v.refPositions,
			MinSize:        v.minSize,
			Append:         v.append,
		}
		faFilter(s)
		records := fasta.Read(v.outputFile)
		expected := fasta.Read(v.expectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error in faFilter. Expected: %s, Out: %s", v.expectedFile, v.outputFile)
		}
		err = os.Remove(v.outputFile)
		exception.PanicOnErr(err)
	}
}
