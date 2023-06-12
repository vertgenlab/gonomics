package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

var FaFormatTests = []struct {
	inputFile        string
	outputFile       string
	expectedFile     string
	lineLength       int
	nameFile         string
	trimName         bool
	toUpper          bool
	revComp          bool
	noGaps           bool
	noGapBed         string
	noGapBedExpected string
	maskInvalid      bool
}{
	{"testdata/faFormatTest.fa", "testdata/faFormatOutput.fa", "testdata/faFormatExpected.fa", 50, "", true, true, false, true, "testdata/test.NoGap.bed", "testdata/expected.NoGap.bed", false},
	{"testdata/faFormatTest.fa", "testdata/faFormatOutput.fa", "testdata/faFormatNamesExpected.fa", 50, "testdata/fastaNames.txt", true, true, false, false, "", "", false},
	{"testdata/revCompTest.fa", "testdata/revCompOutput.fa", "testdata/revCompExpected.fa", 50, "", false, false, true, false, "", "", false},
	{"testdata/revCompTest.fa", "testdata/revCompNamesOutput.fa", "testdata/revCompNamesExpected.fa", 50, "testdata/fastaNames.txt", false, false, true, false, "", "", false},
	{"testdata/maskInput.fa", "testdata/maskOutput.fa", "testdata/maskExpected.fa", 19, "", false, false, false, false, "", "", true},
}

func TestFaFormat(t *testing.T) {
	var err error
	for _, v := range FaFormatTests {
		s := Settings{
			InFile:      v.inputFile,
			OutFile:     v.outputFile,
			LineLength:  v.lineLength,
			NamesFile:   v.nameFile,
			TrimName:    v.trimName,
			ToUpper:     v.toUpper,
			RevComp:     v.revComp,
			NoGaps:      v.noGaps,
			NoGapBed:    v.noGapBed,
			MaskInvalid: v.maskInvalid,
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
