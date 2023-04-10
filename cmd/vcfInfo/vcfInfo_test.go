package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var VcfInfoTests = []struct {
	InFile                 string
	TypesOutFile           string
	DivergenceOutFile      string
	TypesExpectedFile      string
	DivergenceExpectedFile string
}{
	{"testdata/test.vcf", "testdata/tmpTypes.txt", "", "testdata/expectedTypes.txt", ""},
	{"testdata/test.vcf", "", "testdata/tmpDiverge.txt", "", "testdata/expectedDiverge.txt"},
}

func TestVcfInfo(t *testing.T) {
	var err error
	for _, v := range VcfInfoTests {
		handleInputs(v.InFile, v.TypesOutFile, v.DivergenceOutFile, "", false, 0, "", 0)
		if v.TypesOutFile != "" && !fileio.AreEqual(v.TypesExpectedFile, v.TypesOutFile) {
			t.Errorf("Error in vcfInfo. Output file did not match expected.")
		} else if v.TypesOutFile != "" {
			err = os.Remove(v.TypesOutFile)
			exception.PanicOnErr(err)
		}
		if v.DivergenceOutFile != "" && !fileio.AreEqual(v.DivergenceExpectedFile, v.DivergenceOutFile) {
			t.Errorf("Error in vcfInfo. Output file did not match expected.")
		} else if v.DivergenceOutFile != "" {
			err = os.Remove(v.DivergenceOutFile)
			exception.PanicOnErr(err)
		}
	}
}

func TestVcfContext(t *testing.T) {
	fa := "testdata/test.fasta"
	v := "testdata/testContext.vcf"
	expectedMerge := "testdata/expectedMergeComplements.txt"
	expectedInclude := "testdata/expectedIncludeComplements.txt"
	outMerge := "testdata/outMergeComplements.txt"
	outInclude := "testdata/outIncludeComplements.txt"
	var err error

	handleInputs(v, "", "", outMerge, false, 1, fa, 0)
	handleInputs(v, "", "", outInclude, true, 1, fa, 0)

	if !fileio.AreEqual(expectedMerge, outMerge) || !fileio.AreEqual(expectedInclude, outInclude) {
		t.Error("ERROR: problem with vcfContext")
		return
	}

	err = os.Remove(outMerge)
	exception.PanicOnErr(err)
	err = os.Remove(outInclude)
	exception.PanicOnErr(err)
}
