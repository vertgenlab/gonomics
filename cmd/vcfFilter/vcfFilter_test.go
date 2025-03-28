package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"strings"
	"testing"
)

var VcfFilterTests = []struct {
	inputFile                      string
	tmpOutFile                     string
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
	minDaf                         float64
	maxDaf                         float64
	minDist                        int
	setSeed                        int64
}{
	{"testdata/test.vcf", "testdata/tmp.Out.vcf", "testdata/expectedOut.vcf", "testdata/test.group", "chr3", 10, 1000, 0, "", "", true, true, true, false, false, false, false, false, false, false, false, "", 1, 0, 1, 0, 10},
	{"testdata/test_removeNoAncestor.vcf", "testdata/tmp.removeNoAncestor.vcf", "testdata/expected_removeNoAncestor.vcf", "", "", 0, 100, 0, "", "", false, false, false, true, false, false, false, false, false, false, false, "", 1, 0, 1, 0, 10},
	{"testdata/test_onlyPolarizable.vcf", "testdata/tmp.OnlyPolarizable.vcf", "testdata/expected_onlyPolarizable.vcf", "", "", 0, 100, 0, "", "", false, false, false, false, true, false, false, false, false, false, false, "", 1, 0, 1, 0, 10},
	{"testdata/test_weakToStrong.vcf", "testdata/tmp.weakToStrong.vcf", "testdata/expected_noWeakToStrongOrStrongToWeak.vcf", "", "", 0, 100, 0, "", "", false, false, false, false, false, false, true, false, false, false, false, "", 1, 0, 1, 0, 10},
	{"testdata/test_weakToStrong.vcf", "tmp.refWeakAltStrong.vcf", "testdata/expected_refWeakAltStrongOnly.vcf", "", "", 0, 100, 0, "", "", false, false, false, false, false, false, false, true, false, false, false, "", 1, 0, 1, 0, 10},
	{"testdata/test_id.vcf", "testdata/tmp.id.vcf", "testdata/expected_id.vcf", "testdata/test.group", "chr3", 10, 1000, 0, "", "", false, true, true, false, false, false, false, false, false, false, false, "TestingId", 1, 0, 1, 0, 10},
	{"testdata/test.vcf", "testdata/tmp.subset.vcf", "testdata/expectedSubSet.vcf", "", "chr3", 0, numbers.MaxInt, 0, "", "", false, false, false, false, false, false, false, false, false, false, false, "", 0.5, 0, 1, 0, 20},
	{"testdata/testDuplicatePos.vcf", "testdata/tmp.duplicatePos.vcf", "testdata/expectedDuplicatePos.vcf", "", "", 0, numbers.MaxInt, 0, "", "", true, false, false, false, false, false, false, false, false, false, false, "", 1, 0, 1, 0, 10},
	{"testdata/testDaf.vcf", "testdata/tmp.minMaxDaf.vcf", "testdata/expectedMinMaxDaf.vcf", "", "", 0, numbers.MaxInt, 0, "", "", false, false, false, false, false, false, false, false, false, false, false, "", 1, 0.25, 0.75, 0, 10},
	{"testdata/testMinDist.vcf", "testdata/tmp.minDist.vcf", "testdata/expectedMinDist.vcf", "", "", 0, numbers.MaxInt, 0, "", "", false, false, false, false, false, false, false, false, false, false, false, "", 1, 0, 1, 100, 10},
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
			minDaf:                         v.minDaf,
			maxDaf:                         v.maxDaf,
			subSet:                         v.subSet,
			minDist:                        v.minDist,
		}

		vcfFilter(v.inputFile, v.tmpOutFile, c, v.groupFile, false, false, v.setSeed)
		records, _ := vcf.Read(v.tmpOutFile)
		expected, _ := vcf.Read(v.expectedOutputFile)
		if !vcf.AllEqual(records, expected) {
			t.Errorf("Error in vcfFilter. Tmp: %v. Expected: %v.", v.tmpOutFile, v.expectedOutputFile)
		} else {
			err = os.Remove(v.tmpOutFile)
		}
		exception.PanicOnErr(err)
	}
}
