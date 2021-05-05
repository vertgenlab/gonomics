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
}{
	{"testdata/trimNameTest.fa", "testdata/trimNameOutput.fa", "testdata/trimNameExpected.fa", 50, true},
}

func TestFaFormat(t *testing.T) {
	var err error
	for _, v := range FaFormatTests {
		faFormat(v.inputFile, v.outputFile, v.lineLength, v.trimName)
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
