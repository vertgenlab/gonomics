package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"strings"
	"testing"
)

var VcfFilterTests = []struct {
	inputFile                      string
	expectedOutputFile             string
	groupFile                      string
	chrom                          string
	minPos                         int
	maxPos                         int
	minQual                        float64
	ref                            string
	alt                            string
	biAllelicOnly                  bool
	substitutionsOnly              bool
	segregatingSitesOnly           bool
	removeNoAncestor               bool
	onlyPolarizableAncestors       bool
	weakToStrongOrStrongToWeakOnly bool
	noWeakToStrongOrStrongToWeak   bool
	refWeakAltStrongOnly           bool
	refStrongAltWeakOnly           bool
	notRefStrongAltWeak            bool
	notRefWeakAltStrong            bool
	id                             string
	subSet                         float64
}{
	{"testdata/test.vcf", "testdata/expectedOut.vcf", "testdata/test.group", "chr3", 10, 1000, 0, "", "", true, true, true, false, false, false, false, false, false, false, false, "", 1},
	{"testdata/test_removeNoAncestor.vcf", "testdata/expected_removeNoAncestor.vcf", "", "", 0, 100, 0, "", "", false, false, false, true, false, false, false, false, false, false, false, "", 1},
	{"testdata/test_onlyPolarizable.vcf", "testdata/expected_onlyPolarizable.vcf", "", "", 0, 100, 0, "", "", false, false, false, false, true, false, false, false, false, false, false, "", 1},
	{"testdata/test_weakToStrong.vcf", "testdata/expected_noWeakToStrongOrStrongToWeak.vcf", "", "", 0, 100, 0, "", "", false, false, false, false, false, false, true, false, false, false, false, "", 1},
	{"testdata/test_weakToStrong.vcf", "testdata/expected_refWeakAltStrongOnly.vcf", "", "", 0, 100, 0, "", "", false, false, false, false, false, false, false, true, false, false, false, "", 1},
	{"testdata/test_id.vcf", "testdata/expected_id.vcf", "testdata/test.group", "chr3", 10, 1000, 0, "", "", true, true, true, false, false, false, false, false, false, false, false, "TestingId", 1},
	{"testdata/test.vcf", "testdata/expectedSubSet.vcf", "", "chr3", 0, numbers.MaxInt, 0, "", "", false, false, false, false, false, false, false, false, false, false, false, "", 0.5},
}

func TestVcfFilter(t *testing.T) {
	var err error
	for _, v := range VcfFilterTests {

		var altSlice []string
		if v.alt != "" {
			altSlice = strings.Split(v.alt, ",")
		}
		c := criteria{
			chrom:                          v.chrom,
			groupFile:                      v.groupFile,
			minPos:                         v.minPos,
			maxPos:                         v.maxPos,
			minQual:                        v.minQual,
			ref:                            v.ref,
			alt:                            altSlice,
			biAllelicOnly:                  v.biAllelicOnly,
			substitutionsOnly:              v.substitutionsOnly,
			segregatingSitesOnly:           v.segregatingSitesOnly,
			removeNoAncestor:               v.removeNoAncestor,
			onlyPolarizableAncestors:       v.onlyPolarizableAncestors,
			weakToStrongOrStrongToWeakOnly: v.weakToStrongOrStrongToWeakOnly,
			noWeakToStrongOrStrongToWeak:   v.noWeakToStrongOrStrongToWeak,
			refWeakAltStrongOnly:           v.refWeakAltStrongOnly,
			refStrongAltWeakOnly:           v.refStrongAltWeakOnly,
			notRefStrongAltWeak:            v.notRefStrongAltWeak,
			notRefWeakAltStrong:            v.notRefWeakAltStrong,
			id:                             v.id,
			subSet:                         v.subSet,
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
