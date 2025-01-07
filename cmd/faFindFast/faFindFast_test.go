package main

import (
	"math"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var FaFindFastTests = []struct {
	InFile          string
	OutFile         string
	ExpectedFile    string
	FirstQueryName  string
	SecondQueryName string
	WindowSize      int
	RefChromName    string
	RemoveN         bool
	DivergenceRate  float64
	LongOutput      bool
	OutputAlnPos    bool
	GcContent       bool
}{
	{InFile: "testdata/test_indel.fa", //also test for extra species here
		OutFile:         "testdata/tmp.out.bed",
		ExpectedFile:    "testdata/expected.bed",
		FirstQueryName:  "Human",
		SecondQueryName: "Chimp",
		WindowSize:      10,
		RefChromName:    "chr1",
		RemoveN:         false,
		DivergenceRate:  math.MaxFloat64,
		LongOutput:      false,
		OutputAlnPos:    false,
		GcContent:       false},
	{InFile: "testdata/test_indel.fa",
		OutFile:         "testdata/tmp.noN.bed",
		ExpectedFile:    "testdata/expected.noN.bed",
		FirstQueryName:  "", //also test for FirstQueryName and SecondQueryName defaults here
		SecondQueryName: "",
		WindowSize:      10,
		RefChromName:    "chr1",
		RemoveN:         true,
		DivergenceRate:  math.MaxFloat64,
		LongOutput:      false,
		OutputAlnPos:    false,
		GcContent:       false},
	{InFile: "testdata/test.endDoubleGaps.fa",
		OutFile:         "testdata/tmp.doubleGaps.bed",
		ExpectedFile:    "testdata/expected.bed",
		FirstQueryName:  "Human",
		SecondQueryName: "Gorilla", //also test for a different QueryName here
		WindowSize:      10,
		RefChromName:    "chr1",
		RemoveN:         false,
		DivergenceRate:  math.MaxFloat64,
		LongOutput:      false,
		OutputAlnPos:    false,
		GcContent:       false},
	{InFile: "testdata/test.endGapsQuery.fa",
		OutFile:         "testdata/tmp.queryGaps.bed",
		ExpectedFile:    "testdata/expected.endGapsQuery.bed",
		FirstQueryName:  "Human",
		SecondQueryName: "Chimp",
		WindowSize:      10,
		RefChromName:    "chr1",
		RemoveN:         false,
		DivergenceRate:  math.MaxFloat64,
		LongOutput:      false,
		OutputAlnPos:    false,
		GcContent:       false},
	{InFile: "testdata/test.endGapsRef.fa",
		OutFile:         "testdata/tmp.refGaps.bed",
		ExpectedFile:    "testdata/expected.endGapsRef.bed",
		FirstQueryName:  "Human", //also test for finding a FirstQueryName that is not the 1st sequence
		SecondQueryName: "Chimp",
		WindowSize:      10,
		RefChromName:    "chr1",
		RemoveN:         false,
		DivergenceRate:  math.MaxFloat64,
		LongOutput:      false,
		OutputAlnPos:    false,
		GcContent:       false},
	{InFile: "testdata/test.endGapsRef.fa",
		OutFile:         "testdata/tmp.longOutput.bed",
		ExpectedFile:    "testdata/expected.longOutput.bed",
		FirstQueryName:  "Human",
		SecondQueryName: "Chimp",
		WindowSize:      10,
		RefChromName:    "chr1",
		RemoveN:         false,
		DivergenceRate:  0.01,
		LongOutput:      true,
		OutputAlnPos:    false,
		GcContent:       false},
	{InFile: "testdata/test_indel_3seq.fa",
		OutFile:         "testdata/tmp.out.3seq.bed",
		ExpectedFile:    "testdata/expected_3seq.bed",
		FirstQueryName:  "HumanFirstQuery",
		SecondQueryName: "ChimpSecondQuery",
		WindowSize:      10,
		RefChromName:    "chr1",
		RemoveN:         false,
		DivergenceRate:  math.MaxFloat64,
		LongOutput:      false,
		OutputAlnPos:    false,
		GcContent:       false},
	{InFile: "testdata/test_indel_3seq.fa",
		OutFile:         "testdata/tmp.AlnPos.bed",
		ExpectedFile:    "testdata/expected_alnPos.bed",
		FirstQueryName:  "HumanFirstQuery",
		SecondQueryName: "ChimpSecondQuery",
		WindowSize:      10,
		RefChromName:    "chr1",
		RemoveN:         false,
		DivergenceRate:  math.MaxFloat64,
		LongOutput:      false,
		OutputAlnPos:    true, //test for outputAlnPos
		GcContent:       false},
	{InFile: "testdata/test_indel_3seq.fa",
		OutFile:         "testdata/tmp.out.longAndAlnPos.bed",
		ExpectedFile:    "testdata/expected_longAndAlnPos.bed",
		FirstQueryName:  "HumanFirstQuery",
		SecondQueryName: "ChimpSecondQuery",
		WindowSize:      10,
		RefChromName:    "chr1",
		RemoveN:         false,
		DivergenceRate:  math.MaxFloat64,
		LongOutput:      true, //test for longOutput AND outputAlnPos
		OutputAlnPos:    true,
		GcContent:       false},
}

func TestFaFindFast(t *testing.T) {
	var err error
	for _, v := range FaFindFastTests {
		s := Settings{
			InFile:          v.InFile,
			OutFile:         v.OutFile,
			FirstQueryName:  v.FirstQueryName,
			SecondQueryName: v.SecondQueryName,
			WindowSize:      v.WindowSize,
			RefChromName:    v.RefChromName,
			RemoveN:         v.RemoveN,
			LongOutput:      v.LongOutput,
			DivergenceRate:  v.DivergenceRate,
			OutputAlnPos:    v.OutputAlnPos,
			GcContent:       v.GcContent,
		}
		faFindFast(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in faFindFast. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
