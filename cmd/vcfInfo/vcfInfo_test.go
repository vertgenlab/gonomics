package main

import (
	"testing"
	"os"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/exception"
)

var VcfInfoTests = []struct {
	InFile string
	OutFile string
	ExpectedFile string
	printNumDivergent bool
}{
	{"testdata/test.vcf", "testdata/tmp.txt", "testdata/expected.txt", false},
	{"testdata/test.vcf", "testdata/tmpDiverge.txt", "testdata/expectedDiverge.txt", true},
}

func TestVcfInfo(t *testing.T) {
	var err error
	for _, v := range VcfInfoTests {
		vcfInfo(v.InFile, v.OutFile, v.printNumDivergent)
		if !fileio.AreEqual(v.ExpectedFile, v.OutFile) {
			t.Errorf("Error in vcfInfo. Output file did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
