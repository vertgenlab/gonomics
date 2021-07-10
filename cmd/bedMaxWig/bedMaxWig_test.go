package main

import (
	"github.com/vertgenlab/gonomics/common"
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
}{
	{"testdata/testBed.bed", "testdata/startOneStepOne.wig", "testdata/fake.chrom.sizes", "testdata/testBMWOutput.bed", "testdata/testBMWExpected.bed"},
	//{"testdata/testBed.bed", "testdata/startNonOneStepNonOne.wig", "testdata/fake.chrom.sizes", "testdata/testBMWOutput.bed", "testdata/testBMWExpectedNonOneStartStep.bed"},
	//TODO: If bedMaxWig is written to handle a step other than 1, come back to ensure there is a test for this.
}

func TestBedMaxWig(t *testing.T) {
	var err error
	for _, v := range bedMaxWigTests {
		bedMaxWig(v.inputBed, v.inputWig, v.inputChromSizes, v.outputFile, false)
		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error in bedMaxWig, the output beds is not as expected")
		}
		err = os.Remove(v.outputFile)
		if err != nil {
			common.ExitIfError(err)
		}
	}
}

var bedMaxWigNormFlagTest = []struct {
	inputBed         string
	inputWig         string
	inputChromSizes  string
	outputFileNorm   string
	expectedFileNorm string
}{
	{"testdata/testBed.bed", "testdata/startOneStepOne.wig", "testdata/fake.chrom.sizes", "testdata/testBMWOutputNorm.bed", "testdata/testBMWExpectedNormFlagStep1.bed"},
}

func TestBedMaxWigNormFlag(t *testing.T) {
	var err error
	for _, v := range bedMaxWigNormFlagTest {
		bedMaxWig(v.inputBed, v.inputWig, v.inputChromSizes, v.outputFileNorm, true)
		if !fileio.AreEqual(v.outputFileNorm, v.expectedFileNorm) {
			t.Errorf("Error in bedMaxWig, the output bed with norm flag is not as expected")
		}
		err = os.Remove(v.outputFileNorm) //Remove and see what the file looks like. Change expected file to have dots.
		if err != nil {
			common.ExitIfError(err)
		}
	}
}
