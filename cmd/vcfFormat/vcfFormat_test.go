package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var VcfFormatTests = []struct {
	InFile        string
	OutFile       string
	ExpectedFile  string
	EnsemblToUCSC bool
	UCSCToEnsembl bool
	FixVcfRecords bool
	Ref           string
	ClearInfo     bool
}{
	{"testdata/test.UCSC.vcf", "testdata/tmp.UCSCtoEnsembl.vcf", "testdata/test.Ensembl.vcf", false, true, false, "", false},
	{"testdata/test.Ensembl.vcf", "testdata/tmp.EnsemblToUCSC.vcf", "testdata/test.UCSC.vcf", true, false, false, "", false},
	{"testdata/test.UCSC.vcf", "testdata/tmp.UCSCnoInfo.vcf", "testdata/expected.noInfo.vcf", false, false, false, "", true},
	{"testdata/test.broken.vcf", "testdata/tmp.fixed.vcf", "testdata/expected.fixed.vcf", false, false, true, "testdata/test.fa", false},
}

func TestVcfFormat(t *testing.T) {
	var err error
	for _, v := range VcfFormatTests {
		vcfFormat(v.InFile, v.OutFile, v.EnsemblToUCSC, v.UCSCToEnsembl, v.FixVcfRecords, v.Ref, v.ClearInfo, false)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in VcfFormat. Output does not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
