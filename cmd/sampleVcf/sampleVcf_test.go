package main

import (
	"fmt"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/vcf"
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
		sampleVcf(v.inputFile, "tmp.vcf", v.numVariants, v.numSamples, 0)
		records, recHeader := vcf.Read("tmp.vcf")
		expected, expectedHeader := vcf.Read(v.expectedOutputFile)
		if !vcf.AllEqual(records, expected) {
			fmt.Println(records)
			fmt.Println(expected)
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
