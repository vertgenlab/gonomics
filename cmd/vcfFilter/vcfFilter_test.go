package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"testing"
)

var VcfFilterTests = []struct {
	inputFile            string
	expectedOutputFile   string
	groupFile            string
	chrom                string
	minPos               int
	maxPos               int
	minQual              float64
	ref                  string
	alt                  string
	biAllelicOnly        bool
	substitutionsOnly    bool
	segregatingSitesOnly bool
	removeNoAncestor     bool
	onlyPolarizableAncestors	bool
}{
	{"testdata/test.vcf", "testdata/expectedOut.vcf", "testdata/test.group", "chr3", 10, 1000, 0.0, "", "", true, true, true, false, false},
	{"testdata/test_removeNoAncestor.vcf", "testdata/expected_removeNoAncestor.vcf", "", "", 0, 100, 0.0, "", "", false, false, false, true, false},
	{"testdata/test_onlyPolarizable.vcf", "testdata/expected_onlyPolarizable.vcf", "", "", 0, 100, 0.0, "", "", false, false, false, false, true},
}

func TestVcfFilter(t *testing.T) {
	var err error
	for _, v := range VcfFilterTests {
		vcfFilter(v.inputFile, "tmp.vcf", v.groupFile, v.chrom, v.minPos, v.maxPos, v.ref, v.alt, v.minQual, v.biAllelicOnly, v.substitutionsOnly, v.segregatingSitesOnly, v.removeNoAncestor, v.onlyPolarizableAncestors)
		records := vcf.Read("tmp.vcf")
		expected := vcf.Read(v.expectedOutputFile)
		if !vcf.AllEqual(records, expected) {
			t.Errorf("Error in vcfFilter.")
		}
		err = os.Remove("tmp.vcf")
		if err != nil {
			common.ExitIfError(err)
		}
	}
}
