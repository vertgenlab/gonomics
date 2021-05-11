package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"testing"
	"os"
	"github.com/vertgenlab/gonomics/common"
)

var bedMaxWigTests = []struct {
	inputBed    string
	inputWig	string
	inputChromSizes	string
	outputFile   string
	expectedFile string
}{
	{"testdata/testBed.bed", "testdata/startOneStepOne.wig", "testdata/fake.chrom.sizes", "testdata/testBMWOutput.bed", "testdata/testBMWExpected.bed"},
}


func TestBedMaxWig(t *testing.T) {
	var err error
	for _, v := range bedMaxWigTests {
		bedMaxWig(v.inputBed, v.inputWig, v.inputChromSizes, v.outputFile)
		records := bed.Read(v.outputFile)
		expected := bed.Read(v.expectedFile)
		if !bed.AllAreEqual(records, expected) {
				t.Errorf("Error in bedMaxWig, the output beds is not as expected")
		}
		err = os.Remove(v.outputFile)
		if err != nil {
			common.ExitIfError(err)
		}
	}
}