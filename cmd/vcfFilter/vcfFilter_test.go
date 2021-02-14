package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"testing"
)

var VcfFilterTests = []struct {
	inputFile          string
	expectedOutputFile string
<<<<<<< HEAD
	groupFile string
	chrom string
	minPos int
	maxPos int
	minQual float64
	ref string
	alt string
	biAllelicOnly bool
	substitutionsOnly bool
=======
	groupFile          string
	chrom              string
	minPos             int
	maxPos             int
	minQual            float64
	ref                string
	alt                string
>>>>>>> c86ce339145ec91f4c13224f045bedf538be2a90
}{
	{"testdata/test.vcf", "testdata/expectedOut.vcf", "testdata/test.group", "chr3", 10, 1000, 0.0, "", "", true, true},
}

func TestVcfFilter(t *testing.T) {
	var err error
	for _, v := range VcfFilterTests {
		vcfFilter(v.inputFile, "tmp.vcf", v.groupFile, v.chrom, v.minPos, v.maxPos, v.ref, v.alt, v.minQual, v.biAllelicOnly, v.substitutionsOnly)
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
