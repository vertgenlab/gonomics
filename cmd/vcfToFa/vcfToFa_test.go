package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var VcfToFaTests = []struct {
	inputVcfFile       string
	inputFaFile        string
	actualOutputFile   string
	expectedOutputFile string
}{
	{"testdata/testInput.vcf", "testdata/testInput.fa", "testdata/actual.fa", "testdata/expected.fa"},
}

func TestVcfToFa(t *testing.T) {
	var err error
	for _, v := range VcfToFaTests {
		vcfToFa(v.inputVcfFile, v.inputFaFile, v.actualOutputFile, true)
		if !fileio.AreEqual(v.actualOutputFile, v.expectedOutputFile) {
			t.Errorf("VcfToFa output did not match expected output.")
		} else {
			err = os.Remove(v.actualOutputFile)
			common.ExitIfError(err)
		}
	}
}
