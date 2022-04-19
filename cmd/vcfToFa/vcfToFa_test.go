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
	useAlt             bool
	useSamples         bool
	haploid            bool
	multiFaMode        bool
	multiFaChromName   string
}{
	//{"testdata/testInput.vcf", "testdata/testInput.fa", "testdata/actual.fa", "testdata/expected.fa", true, false, ""},
	//{"testdata/testInput.vcf", "testdata/testInput.fa", "testdata/actualNoAlt.fa", "testdata/testInput.fa", false, false, ""},//should just produce the input file.
	{"testdata/testMultiInput.vcf", "testdata/testMultiInput.fa", "testdata/actualMultiNoAlt.fa", "testdata/expectedMultiNoAlt.fa", false, true, true, true, "chr1"}, //should also produce the input file
	{"testdata/testMultiInput.vcf", "testdata/testMultiInput.fa", "testdata/actualMultiAlt.fa", "testdata/expectedMultiAlt.fa", true, false, false, true, "chr1"},
	{"testdata/testMultiInput.vcf", "testdata/testMultiInput.fa", "testdata/actualMultiAltchr2.fa", "testdata/expectedMultiAltchr2.fa", true, false, false, true, "chr2"},
}

func TestVcfToFa(t *testing.T) {
	var err error
	for _, v := range VcfToFaTests {
		if v.multiFaMode {
			vcfToMultiFa(v.inputVcfFile, v.inputFaFile, v.actualOutputFile, v.multiFaChromName, v.useAlt, v.useSamples, v.haploid)
		} else {
			vcfToFa(v.inputVcfFile, v.inputFaFile, v.actualOutputFile, v.useAlt)
		}
		if !fileio.AreEqual(v.actualOutputFile, v.expectedOutputFile) {
			t.Errorf("VcfToFa output did not match expected output.")
		} else {
			err = os.Remove(v.actualOutputFile)
			common.ExitIfError(err)
		}
	}
}
