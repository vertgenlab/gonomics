package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var MultiFaToVcfTests = []struct {
	inputFa               string
	chromName             string
	outputVcf             string
	expectedVcf           string
	substitutionsOnlyFlag bool
	retainNFlag           bool
}{
	{"testdata/inputMulti.fa", "chr2", "testdata/output.vcf", "testdata/expected.vcf", false, false},
	{"testdata/inputMulti.fa", "chr2", "testdata/output.vcf", "testdata/expectedSubOnly.vcf", true, false},
	{"testdata/inputMulti.fa", "chr2", "testdata/output.vcf", "testdata/expectedRetainN.vcf", false, true},
}

func TestMultiFaVcf(t *testing.T) {
	var err error
	for _, v := range MultiFaToVcfTests {
		multiFaToVcf(v.inputFa, v.chromName, v.outputVcf, v.substitutionsOnlyFlag, v.retainNFlag)
		if !fileio.AreEqual(v.expectedVcf, v.outputVcf) {
			t.Errorf("Error in multiFaToVcf.go, expected.vcf != output.vcf")
		} else {
			err = os.Remove(v.outputVcf)
			exception.PanicOnErr(err)
		}
	}
}
