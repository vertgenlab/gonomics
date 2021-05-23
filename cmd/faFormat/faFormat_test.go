package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"os"
	"testing"
)

var FaFormatTests = []struct {
	inputFile    string
	outputFile   string
	expectedFile string
	lineLength   int
	trimName     bool
	toUpper	bool
	revComp bool
	noGaps bool
}{
	{"testdata/faFormatTest.fa", "testdata/faFormatOutput.fa", "testdata/faFormatExpected.fa", 50, true, true, false, true},
	{"testdata/revCompTest.fa", "testdata/revCompOutput.fa", "testdata/revCompExpected.fa", 50, false, false, true, false},
}

func TestFaFormat(t *testing.T) {
	var err error
	for _, v := range FaFormatTests {
		s := Settings{
			InFile: v.inputFile,
			OutFile: v.outputFile,
			LineLength: v.lineLength,
			TrimName: v.trimName,
			ToUpper: v.toUpper,
			RevComp: v.revComp,
			NoGaps: v.noGaps,
		}
		faFormat(s)
		records := fasta.Read(v.outputFile)
		expected := fasta.Read(v.expectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error in faFormat.")
		}
		err = os.Remove(v.outputFile)
		if err != nil {
			common.ExitIfError(err)
		}
	}
}
