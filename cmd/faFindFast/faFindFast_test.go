package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var FaFindFastTests = []struct {
	InFile         string
	OutFile        string
	ExpectedFile   string
	ReferenceName  string
	QueryName      string
	WindowSize     int
	ChromName      string
	RemoveN        bool
	DivergenceRate float64
	LongOutput     bool
}{
	{InFile: "testdata/test_indel.fa",
		OutFile:        "testdata/tmp.out.bed",
		ExpectedFile:   "testdata/expected.bed",
		ReferenceName:  "Human",
		QueryName:      "Chimp",
		WindowSize:     10,
		ChromName:      "chr1",
		RemoveN:        false,
		DivergenceRate: -1,
		LongOutput:     false},
	{InFile: "testdata/test_indel.fa",
		OutFile:        "testdata/tmp.noN.bed",
		ExpectedFile:   "testdata/expected.noN.bed",
		ReferenceName:  "", //also test for ReferenceName and QueryName defaults here
		QueryName:      "",
		WindowSize:     10,
		ChromName:      "chr1",
		RemoveN:        true,
		DivergenceRate: -1,
		LongOutput:     false},
	{InFile: "testdata/test.endDoubleGaps.fa",
		OutFile:        "testdata/tmp.doubleGaps.bed",
		ExpectedFile:   "testdata/expected.bed",
		ReferenceName:  "Human",
		QueryName:      "Gorilla", //also test for a different QueryName here
		WindowSize:     10,
		ChromName:      "chr1",
		RemoveN:        false,
		DivergenceRate: -1,
		LongOutput:     false},
	{InFile: "testdata/test.endGapsQuery.fa",
		OutFile:        "testdata/tmp.queryGaps.bed",
		ExpectedFile:   "testdata/expected.endGapsQuery.bed",
		ReferenceName:  "Human",
		QueryName:      "Chimp",
		WindowSize:     10,
		ChromName:      "chr1",
		RemoveN:        false,
		DivergenceRate: -1,
		LongOutput:     false},
	{InFile: "testdata/test.endGapsRef.fa",
		OutFile:        "testdata/tmp.refGaps.bed",
		ExpectedFile:   "testdata/expected.endGapsRef.bed",
		ReferenceName:  "Human", //also test for finding a ReferenceName which is not the 1st sequence
		QueryName:      "Chimp",
		WindowSize:     10,
		ChromName:      "chr1",
		RemoveN:        false,
		DivergenceRate: -1,
		LongOutput:     false},
	{InFile: "testdata/test.endGapsRef.fa",
		OutFile:        "testdata/tmp.longOutput.bed",
		ExpectedFile:   "testdata/expected.longOutput.bed",
		ReferenceName:  "Human",
		QueryName:      "Chimp",
		WindowSize:     10,
		ChromName:      "chr1",
		RemoveN:        false,
		DivergenceRate: 0.01,
		LongOutput:     true},
}

func TestFaFindFast(t *testing.T) {
	var err error
	for _, v := range FaFindFastTests {
		faFindFast(v.InFile, v.OutFile, v.ReferenceName, v.QueryName, v.WindowSize, v.ChromName, v.RemoveN, v.LongOutput, v.DivergenceRate)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in faFindFast. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
