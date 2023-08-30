package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

// TODO: fix according to new faFindFast.go

var FaFindFastTests = []struct {
	InFile          string
	OutFile         string
	ExpectedFile    string
	ReferenceName   string
	QueryName       string
	PosRefName      string
	WindowSize      int
	posRefChromName string
	RemoveN         bool
	DivergenceRate  float64
	LongOutput      bool
}{
	{InFile: "testdata/test_indel.fa",
		OutFile:         "testdata/tmp.out.bed",
		ExpectedFile:    "testdata/expected.bed",
		ReferenceName:   "Human", //also test for extra species here
		QueryName:       "Chimp",
		PosRefName:      "Human",
		WindowSize:      10,
		posRefChromName: "chr1",
		RemoveN:         false,
		DivergenceRate:  -1,
		LongOutput:      false},
	{InFile: "testdata/test_indel.fa",
		OutFile:         "testdata/tmp.noN.bed",
		ExpectedFile:    "testdata/expected.noN.bed",
		ReferenceName:   "", //also test for ReferenceName, QueryName and PosRefName defaults here
		QueryName:       "",
		PosRefName:      "",
		WindowSize:      10,
		posRefChromName: "chr1",
		RemoveN:         true,
		DivergenceRate:  -1,
		LongOutput:      false},
	{InFile: "testdata/test.endDoubleGaps.fa",
		OutFile:         "testdata/tmp.doubleGaps.bed",
		ExpectedFile:    "testdata/expected.bed",
		ReferenceName:   "Human",
		QueryName:       "Gorilla", //also test for a different QueryName here
		PosRefName:      "",
		WindowSize:      10,
		posRefChromName: "chr1",
		RemoveN:         false,
		DivergenceRate:  -1,
		LongOutput:      false},
	{InFile: "testdata/test.endGapsQuery.fa",
		OutFile:         "testdata/tmp.queryGaps.bed",
		ExpectedFile:    "testdata/expected.endGapsQuery.bed",
		ReferenceName:   "Human",
		QueryName:       "Chimp",
		PosRefName:      "",
		WindowSize:      10,
		posRefChromName: "chr1",
		RemoveN:         false,
		DivergenceRate:  -1,
		LongOutput:      false},
	{InFile: "testdata/test.endGapsRef.fa",
		OutFile:         "testdata/tmp.refGaps.bed",
		ExpectedFile:    "testdata/expected.endGapsRef.bed",
		ReferenceName:   "Human", //also test for finding a ReferenceName that is not the 1st sequence
		QueryName:       "Chimp",
		PosRefName:      "HumanPosRef", //also test for finding a PosRefName that is not the same as ReferenceName
		WindowSize:      10,
		posRefChromName: "chr1",
		RemoveN:         false,
		DivergenceRate:  -1,
		LongOutput:      false},
	{InFile: "testdata/test.endGapsRef.fa",
		OutFile:         "testdata/tmp.longOutput.bed",
		ExpectedFile:    "testdata/expected.longOutput.bed",
		ReferenceName:   "Human",
		QueryName:       "Chimp",
		PosRefName:      "",
		WindowSize:      10,
		posRefChromName: "chr1",
		RemoveN:         false,
		DivergenceRate:  0.01,
		LongOutput:      true},
}

func TestFaFindFast(t *testing.T) {
	var err error
	for _, v := range FaFindFastTests {
		faFindFast(v.InFile, v.OutFile, v.ReferenceName, v.QueryName, v.PosRefName, v.WindowSize, v.posRefChromName, v.RemoveN, v.LongOutput, v.DivergenceRate)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in faFindFast. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
