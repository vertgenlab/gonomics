package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var SubsequenceSwapTestCases = []struct {
	InFile         string
	BackgroundName string
	ForegroundName string
	SwapRegions    string
	ChromName      string
	OutFile        string
	ExpectedFile   string
}{
	{
		InFile:         "testdata/test.fa",
		BackgroundName: "Seq4",
		ForegroundName: "Seq3",
		SwapRegions:    "testdata/swapRegionsTest.bed",
		ChromName:      "chr1",
		OutFile:        "testdata/outputFile.fa",
		ExpectedFile:   "testdata/expectedSwap.fa",
	},
	{
		InFile:         "testdata/testWithIndels.fa",
		BackgroundName: "hg38",
		ForegroundName: "hca",
		SwapRegions:    "testdata/swapWithIndels.bed",
		ChromName:      "chr1",
		OutFile:        "testdata/test.OutputWithIndels.fa",
		ExpectedFile:   "testdata/expected.SwapWithIndel.fa",
	},
	{
		InFile:         "testdata/test.fa",
		BackgroundName: "Seq1",
		ForegroundName: "Seq2",
		SwapRegions:    "testdata/swapWithChrom.bed",
		ChromName:      "chr1",
		OutFile:        "testdata/test.OutputWithChrom.fa",
		ExpectedFile:   "testdata/expected.SwapWithChrom.fa",
	},
}

func TestSubsequenceSwap(t *testing.T) {
	var err error
	for _, tc := range SubsequenceSwapTestCases {
		//loop through all test cases
		multiFaSubsequenceSwap(tc.InFile, tc.SwapRegions, tc.BackgroundName, tc.ForegroundName, tc.ChromName, tc.OutFile)

		if !fileio.AreEqual(tc.OutFile, tc.ExpectedFile) {
			t.Errorf("Error: Output was not as expected. \n")
		} else {
			err = os.Remove(tc.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
