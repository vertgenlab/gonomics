package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var FaFormatTests = []struct {
	inputFile    string
	outputFile   string
	expectedFile string
	lineLength   int
	trimName     bool
	toUpper      bool
	revComp      bool
	noGaps       bool
	noGapBed	string
	noGapBedExpected string
}{
	{"testdata/faFormatTest.fa", "testdata/faFormatOutput.fa", "testdata/faFormatExpected.fa", 50, true, true, false, true, "testdata/test.NoGap.bed", "testdata/expected.NoGap.bed"},
	{"testdata/revCompTest.fa", "testdata/revCompOutput.fa", "testdata/revCompExpected.fa", 50, false, false, true, false, "", ""},
}

func TestFaFormat(t *testing.T) {
	var err error
	for _, v := range FaFormatTests {
		s := Settings{
			InFile:     v.inputFile,
			OutFile:    v.outputFile,
			LineLength: v.lineLength,
			TrimName:   v.trimName,
			ToUpper:    v.toUpper,
			RevComp:    v.revComp,
			NoGaps:     v.noGaps,
			NoGapBed:	v.noGapBed,
		}
		faFormat(s)
		records := fasta.Read(v.outputFile)
		expected := fasta.Read(v.expectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error in faFormat.")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
		if v.noGapBed != "" {
			if !fileio.AreEqual(v.noGapBed, v.noGapBedExpected) {
				t.Errorf("Error in faFormat, noGapBed did not match expected.")
			} else {
				err = os.Remove(v.noGapBed)
				exception.PanicOnErr(err)
			}
		}
	}
}
