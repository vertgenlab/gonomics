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
}{
	{"testdata/testBed.bed", "testdata/startOneStepOne.wig", "testdata/fake.chrom.sizes", "testdata/testBMWOutput.bed", "testdata/testBMWExpected.bed"},
	//{"testdata/testBed.bed", "testdata/startNonOneStepNonOne.wig", "testdata/fake.chrom.sizes", "testdata/testBMWOutput.bed", "testdata/testBMWExpectedNonOneStartStep.bed"},
	//TODO: If bedMaxWig is written to handle a step other than 1, come back to ensure there is a test for this.
}

func TestBedMaxWig(t *testing.T) {
	var err error
	for _, v := range bedMaxWigTests {
		bedValueWig(v.inputBed, v.inputWig, v.inputChromSizes, v.outputFile, false, false)
		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error in bedMaxWig, the output beds is not as expected")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
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
		bedValueWig(v.inputBed, v.inputWig, v.inputChromSizes, v.outputFileNorm, true, false)
		if !fileio.AreEqual(v.outputFileNorm, v.expectedFileNorm) {
			t.Errorf("Error in bedMaxWig, the output bed with norm flag is not as expected")
		} else {
			err = os.Remove(v.outputFileNorm) //Remove and see what the file looks like. Change expected file to have dots.
			exception.PanicOnErr(err)
		}
	}
}

var bedMinWigTest = []struct {
	inputBed        string
	inputWig        string
	inputChromSizes string
	outputFile      string
	expectedFile    string
}{
	{"testdata/testBed.bed", "testdata/startOneStepOne.wig", "testdata/fake.chrom.sizes", "testdata/testMinOutput.bed", "testdata/testMinExpected.bed"},
}

func TestBedMinWig(t *testing.T) {
	var err error
	for _, v := range bedMinWigTest {
		bedValueWig(v.inputBed, v.inputWig, v.inputChromSizes, v.outputFile, false, true)
		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error in bedMinWig, the output bed is not as expected.")
		} else {
			err = os.Remove(v.outputFile)
			exception.PanicOnErr(err)
		}
	}
}
