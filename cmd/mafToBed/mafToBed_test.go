package main

import (
	"testing"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
)


var MafToBedTests = []struct {
	inputFile string
	outputFile string
	expectedFile string
	reference string
}{
	{"testdata/chr22.test.maf", "testdata/tmp.bed", "testdata/expected.bed", "hg38"},
}

func TestMafToBed(t *testing.T) {
	var err error
	for _, v := range MafToBedTests {
		mafToBed(v.inputFile, v.outputFile, v.reference)
		if !fileio.AreEqual(v.outputFile, v.expectedFile) {
			t.Errorf("Error in mafToBed.")
		} else {
			err = os.Remove(v.outputFile)
			if err != nil {
				common.ExitIfError(err)
			}
		}
	}
}