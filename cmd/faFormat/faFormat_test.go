package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

var FaFormatTests = []struct {
	InputFile               string
	OutputFile              string
	ExpectedFile            string
	LineLength              int
	NameFile                string
	TrimName                bool
	ToUpper                 bool
	ToLower                 string
	RevComp                 bool
	NoGaps                  bool
	NoGapBed                string
	NoGapBedExpected        string
	MaskInvalid             bool
	MultiFaNoGapBed         string
	QuerySeqName            string
	ChromName               string
	ExpectedMultiFaNoGapBed string
	Rename                  string
}{
	{InputFile: "testdata/faFormatTest.fa",
		OutputFile:       "testdata/faFormatOutput.fa",
		ExpectedFile:     "testdata/faFormatExpected.fa",
		LineLength:       50,
		NameFile:         "",
		TrimName:         true,
		ToUpper:          true,
		ToLower:          "",
		RevComp:          false,
		NoGaps:           true,
		NoGapBed:         "testdata/test.NoGap.bed",
		NoGapBedExpected: "testdata/expected.NoGap.bed",
		MaskInvalid:      false},
	{InputFile: "testdata/faFormatTest.fa",
		OutputFile:       "testdata/faFormatOutput.fa",
		ExpectedFile:     "testdata/faFormatNamesExpected.fa",
		LineLength:       50,
		NameFile:         "testdata/fastaNames.txt",
		TrimName:         true,
		ToUpper:          true,
		ToLower:          "",
		RevComp:          false,
		NoGaps:           false,
		NoGapBed:         "",
		NoGapBedExpected: "",
		MaskInvalid:      false},
	{InputFile: "testdata/revCompTest.fa",
		OutputFile:       "testdata/revCompOutput.fa",
		ExpectedFile:     "testdata/revCompExpected.fa",
		LineLength:       50,
		NameFile:         "",
		TrimName:         false,
		ToUpper:          false,
		ToLower:          "",
		RevComp:          true,
		NoGaps:           false,
		NoGapBed:         "",
		NoGapBedExpected: "",
		MaskInvalid:      false},
	{InputFile: "testdata/revCompTest.fa",
		OutputFile:       "testdata/revCompNamesOutput.fa",
		ExpectedFile:     "testdata/revCompNamesExpected.fa",
		LineLength:       50,
		NameFile:         "testdata/fastaNames.txt",
		TrimName:         false,
		ToUpper:          false,
		ToLower:          "",
		RevComp:          true,
		NoGaps:           false,
		NoGapBed:         "",
		NoGapBedExpected: "",
		MaskInvalid:      false},
	{InputFile: "testdata/maskInput.fa",
		OutputFile:       "testdata/maskOutput.fa",
		ExpectedFile:     "testdata/maskExpected.fa",
		LineLength:       19,
		NameFile:         "",
		TrimName:         false,
		ToUpper:          false,
		ToLower:          "",
		RevComp:          false,
		NoGaps:           false,
		NoGapBed:         "",
		NoGapBedExpected: "",
		MaskInvalid:      true},
	{InputFile: "testdata/multiFaGaps.fa",
		OutputFile:              "testdata/out.multiFaGaps.fa",
		ExpectedFile:            "testdata/expected.multiFaGaps.fa",
		LineLength:              50,
		NameFile:                "",
		TrimName:                false,
		ToUpper:                 false,
		ToLower:                 "",
		RevComp:                 false,
		NoGaps:                  false,
		NoGapBed:                "",
		NoGapBedExpected:        "",
		MaskInvalid:             false,
		MultiFaNoGapBed:         "testdata/out.multiFaNoGap.bed",
		QuerySeqName:            "hca",
		ChromName:               "chr1",
		ExpectedMultiFaNoGapBed: "testdata/expected.multiFaNoGap.bed",
	},
	{InputFile: "testdata/faFormatTest.fa",
		OutputFile:       "testdata/out.Rename.fa",
		ExpectedFile:     "testdata/expected.Rename.fa",
		LineLength:       50,
		NameFile:         "",
		TrimName:         false,
		ToUpper:          false,
		ToLower:          "",
		RevComp:          false,
		NoGaps:           false,
		NoGapBed:         "",
		NoGapBedExpected: "",
		MaskInvalid:      false,
		Rename:           "NoGapTest,RenamedField",
	},
	{InputFile: "testdata/toLower.fa",
		OutputFile:       "testdata/out.toLower.fa",
		ExpectedFile:     "testdata/expected.toLower.fa",
		LineLength:       50,
		NameFile:         "",
		TrimName:         false,
		ToUpper:          false,
		ToLower:          "testdata/toLower.bed",
		RevComp:          false,
		NoGaps:           false,
		NoGapBed:         "",
		NoGapBedExpected: "",
		MaskInvalid:      false,
		Rename:           "",
	},
}

func TestFaFormat(t *testing.T) {
	var err error
	for _, v := range FaFormatTests {
		s := Settings{
			InFile:          v.InputFile,
			OutFile:         v.OutputFile,
			LineLength:      v.LineLength,
			NamesFile:       v.NameFile,
			TrimName:        v.TrimName,
			ToUpper:         v.ToUpper,
			ToLower:         v.ToLower,
			RevComp:         v.RevComp,
			NoGaps:          v.NoGaps,
			NoGapBed:        v.NoGapBed,
			MaskInvalid:     v.MaskInvalid,
			MultiFaNoGapBed: v.MultiFaNoGapBed,
			QuerySeqName:    v.QuerySeqName,
			ChromName:       v.ChromName,
			Rename:          v.Rename,
		}
		faFormat(s)
		records := fasta.Read(v.OutputFile)
		expected := fasta.Read(v.ExpectedFile)
		if !fasta.AllAreEqual(records, expected) {
			t.Errorf("Error: in faFormat.")
		} else {
			err = os.Remove(v.OutputFile)
			exception.PanicOnErr(err)
		}
		if v.NoGapBed != "" {
			if !fileio.AreEqual(v.NoGapBed, v.NoGapBedExpected) {
				t.Errorf("Error: in faFormat, noGapBed did not match expected.")
			} else {
				err = os.Remove(v.NoGapBed)
				exception.PanicOnErr(err)
			}
		}
		if v.MultiFaNoGapBed != "" {
			if !fileio.AreEqual(v.MultiFaNoGapBed, v.ExpectedMultiFaNoGapBed) {
				t.Errorf("Error: in faFormat, MultiFaNoGapBed did not match expected.")
			} else {
				err = os.Remove(v.MultiFaNoGapBed)
				exception.PanicOnErr(err)
			}
		}
	}
}
