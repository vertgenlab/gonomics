package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var MultiFaToVcfTests = []struct {
	inputFa               string
	chromName             string
	outputVcf             string
	expectedVcf           string
	substitutionsOnlyFlag bool
	retainNFlag           bool
	secondQueryName       string
}{
	{"testdata/inputMulti.fa", "chr2", "testdata/output.vcf", "testdata/expected.vcf", false, false, ""},
	{"testdata/inputMulti.fa", "chr2", "testdata/output.vcf", "testdata/expectedSubOnly.vcf", true, false, ""},
	{"testdata/inputMulti.fa", "chr2", "testdata/output.vcf", "testdata/expectedRetainN.vcf", false, true, ""},
	{"testdata/inputStartWithGap.fa", "chr2", "testdata/startGapTmp.vcf", "testdata/expectedStartGap.vcf", false, false, ""},
	{"testdata/inputAltStartWithGap.fa", "chr2", "testdata/startGapAltTmp.vcf", "testdata/expectedAltStartsWithGap.vcf", false, false, ""}, //TODO: This expected file does not contain a deletion at position 1, as our parser cannot handle this case.
	{"testdata/inputMultiSecondQueryName.fa", "chr2", "testdata/secondQueryNameTmp.vcf", "testdata/expected.vcf", false, false, "HCA"},
}

func TestMultiFaVcf(t *testing.T) {
	var err error
	for _, v := range MultiFaToVcfTests {
		multiFaToVcf(v.inputFa, v.chromName, v.outputVcf, v.substitutionsOnlyFlag, v.retainNFlag, v.secondQueryName)
		if !fileio.AreEqual(v.expectedVcf, v.outputVcf) {
			t.Errorf("Error in multiFaToVcf.go, expected.vcf != output.vcf")
		} else {
			err = os.Remove(v.outputVcf)
			exception.PanicOnErr(err)
		}
	}
}
