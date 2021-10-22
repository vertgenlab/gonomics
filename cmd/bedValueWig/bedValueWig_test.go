package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var bedMaxWigTests = []struct {
	inputBed        string
	inputWig        string
	inputChromSizes string
	outputFile      string
	expectedFile    string
	minFlag bool
	averageFlag bool
	normFlag bool
}{
	{"testdata/testBed.bed", "testdata/startOneStepOne.wig", "testdata/fake.chrom.sizes", "testdata/testBMWOutput.bed", "testdata/testBMWExpected.bed", false, false, false},
	{"testdata/testBed.bed", "testdata/startOneStepOne.wig", "testdata/fake.chrom.sizes", "testdata/testBMWOutputNorm.bed", "testdata/testBMWExpectedNormFlagStep1.bed", false, false, true},
	{"testdata/testBed.bed", "testdata/startOneStepOne.wig", "testdata/fake.chrom.sizes", "testdata/testMinOutput.bed", "testdata/testMinExpected.bed", true, false, false},
	{"testdata/testBed.bed", "testdata/startOneStepOne.wig", "testdata/fake.chrom.sizes", "testdata/testAverageOutput.bed", "testdata/testAverageExpected.bed", false, true, false},
}

func TestBedValueWig(t *testing.T) {
	var err error
	for _, v := range bedMaxWigTests {
		s := Settings {
			Infile: v.inputBed,
			WigFile: v.inputWig,
			SizesFile: v.inputChromSizes,
			OutFile: v.outputFile,
			MinFlag: v.minFlag,
			AverageFlag: v.averageFlag,
			NormFlag: v.normFlag,
		}
		bedValueWig(s)
		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error in bedMaxWig, the output beds is not as expected")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
	}
}
