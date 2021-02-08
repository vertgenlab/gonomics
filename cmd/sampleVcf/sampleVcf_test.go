package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"testing"
)

var SampleVcfTests = []struct {
	inputFile          string
	expectedOutputFile string
	numVariants        int
	numSamples         int
}{
	{"testdata/babyTest.vcf", "testdata/babyTest_expected.vcf", 3, 2},
}

func TestSampleVcf(t *testing.T) {
	for _, v := range SampleVcfTests {
		sampleVcf(v.inputFile, "tmp.vcf", v.numVariants, v.numSamples, false, 0)
		records, recHeader := vcf.ReadWithHeader("tmp.vcf")
		expected, expectedHeader := vcf.ReadWithHeader(v.expectedOutputFile)
		if !vcf.AllEqual(records, expected) {
			t.Errorf("Error in sampleVcf.")
		}
		if vcf.CompareHeader(recHeader, expectedHeader) != 0 {
			t.Errorf("Error in sampleVcf headers.")
		}
		err := os.Remove("tmp.vcf")
		if err != nil {
			common.ExitIfError(err)
		}
	}
}
