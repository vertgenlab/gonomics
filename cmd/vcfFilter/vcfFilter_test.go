package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"strings"
	"testing"
)

var VcfFilterTests = []struct {
	inputFile                string
	expectedOutputFile       string
	groupFile                string
	chrom                    string
	minPos                   int
	maxPos                   int
	minQual                  float64
	ref                      string
	alt                      string
	biAllelicOnly            bool
	substitutionsOnly        bool
	segregatingSitesOnly     bool
	removeNoAncestor         bool
	onlyPolarizableAncestors bool
	id											 string //raven's note: added id, and added field in corresponding tests below
}{
	{"testdata/test.vcf", "testdata/expectedOut.vcf", "testdata/test.group", "chr3", 10, 1000, 0, "", "", true, true, true, false, false, ""},
	{"testdata/test_removeNoAncestor.vcf", "testdata/expected_removeNoAncestor.vcf", "", "", 0, 100, 0, "", "", false, false, false, true, false, ""},
	{"testdata/test_onlyPolarizable.vcf", "testdata/expected_onlyPolarizable.vcf", "", "", 0, 100, 0, "", "", false, false, false, false, true, ""},
	{"testdata/test_id.vcf", "testdata/expected_id.vcf", "testdata/test.group", "chr3", 10, 1000, 0, "", "", true, true, true, false, false, "TestingId"}, //raven's note: this is the test for id
}

func TestVcfFilter(t *testing.T) {
	var err error
	for _, v := range VcfFilterTests {

		var altSlice []string
		if v.alt != "" {
			altSlice = strings.Split(v.alt, ",")
		}
		c := criteria{
			chrom:                    v.chrom,
			groupFile:                v.groupFile,
			minPos:                   v.minPos,
			maxPos:                   v.maxPos,
			minQual:                  v.minQual,
			ref:                      v.ref,
			alt:                      altSlice,
			biAllelicOnly:            v.biAllelicOnly,
			substitutionsOnly:        v.substitutionsOnly,
			segregatingSitesOnly:     v.segregatingSitesOnly,
			removeNoAncestor:         v.removeNoAncestor,
			onlyPolarizableAncestors: v.onlyPolarizableAncestors,
			id:												v.id, //raven's note: added id
		}

		vcfFilter(v.inputFile, "tmp.vcf", c, v.groupFile, false, false)
		records, _ := vcf.Read("tmp.vcf")
		expected, _ := vcf.Read(v.expectedOutputFile)
		if !vcf.AllEqual(records, expected) {
			t.Errorf("Error in vcfFilter.")
		}
		err = os.Remove("tmp.vcf")
		if err != nil {
			common.ExitIfError(err)
		}
	}
}
