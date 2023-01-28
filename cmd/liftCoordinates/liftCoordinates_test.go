package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"testing"
)

var LiftTests = []struct {
	inputFile          string
	outFile string
	expectedOutputFile string
	unmappedFile string
	chainFile          string
	minMatch float64
	faFile             string
	verbose            int
	swapAB             bool
	strictBorders bool
}{
	{inputFile: "testdata/input.vcf",
		outFile: "testdata/tmp.vcf",
		expectedOutputFile: "testdata/expected.vcf",
		unmappedFile: "testdata/unmapped.txt",
		chainFile: "testdata/test.chain",
		minMatch: 0.95,
		faFile: "testdata/test.fa",
		verbose: 0,
		swapAB: false,
	strictBorders: false},
	{"testdata/input_swapAB.vcf",
		"testdata/testSwap.vcf",
		"testdata/expected_swapAB.vcf",
		"testdata/swap.unmapped.txt",
		"testdata/test.chain",
		0.95,
		"testdata/test.fa",
		0,
		true,
	false},
}

func TestLift(t *testing.T) {
	var err error
	var s Settings

	for _, v := range LiftTests {
		s = Settings {
			InFile:          v.inputFile,
			OutFile: v.outFile,
			UnmappedFile: v.unmappedFile,
			ChainFile:          v.chainFile,
			FaFile:             v.faFile,
			MinMatch: v.minMatch,
			Verbose:           v.verbose,
			SwapAB:             v.swapAB,
			StrictBorders: v.strictBorders,
		}
		liftCoordinates(s)

		if vcf.IsVcfFile(v.inputFile) {
			liftCoordinates(s)
			records, _ := vcf.Read(v.outFile)
			expected, _ := vcf.Read(v.expectedOutputFile)
			if !vcf.AllEqual(records, expected) {
				t.Errorf("Error in Lift for vcf.")
			} else {
				err = os.Remove(v.outFile)
				exception.PanicOnErr(err)
				err = os.Remove(v.unmappedFile)
				exception.PanicOnErr(err)
			}
		} else {
			liftCoordinates(s)
			records := bed.Read("tmp.bed")
			expected := bed.Read(v.expectedOutputFile)
			if !bed.AllAreEqual(records, expected) {
				t.Errorf("Error in Lift for bed.")
			} else {
				err = os.Remove("tmp.bed")
				exception.PanicOnErr(err)
				err = os.Remove("tmp.unmapped")
				exception.PanicOnErr(err)
			}
		}
	}
}
