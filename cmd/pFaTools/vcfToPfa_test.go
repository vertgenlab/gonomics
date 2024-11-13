package main

import (
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

var vcfToPfaTests = []struct {
	InFile	string
	RefFile	string
	OutDir	string
	Start	int
	End		int
	Expected string
	Precision float32
}{
	{InFile: "testdata/test_vcfToPfa_input_1.vcf",
		RefFile:     "testdata/test_vcfToPfa_input_1.fa",
		OutDir:    "testdata/test_vcfToPfa_observed_1.pfa",
		Start:     1,
		End:       36,
		Expected:  "testdata/test_vcfToPfa_expected_1.pfa",
		Precision: 1e-3,
	},
	{InFile: "testdata/test_vcfToPfa_input_2.vcf",
		RefFile:     "testdata/test_vcfToPfa_input_1.fa",
		OutDir:    "testdata/test_vcfToPfa_observed_2.pfa",
		Start:     1,
		End:       44,
		Expected:  "testdata/test_vcfToPfa_expected_2.pfa",
		Precision: 1e-3,
	},
}

func TestVcfToPfa(t *testing.T) {
	for _, tc := range vcfToPfaTests {
		observed := []pFasta.PFasta{pFasta.VcfToPfa(tc.InFile, tc.RefFile, tc.Start, tc.End)}
		pFasta.Write(tc.OutDir, observed)
		expected := pFasta.Read(tc.Expected)
		if !pFasta.AllAreEqual(observed, expected, tc.Precision) {
			t.Errorf("Error: in pFasta. Sample valid input test not as expected.\n")
		} else {
			fileio.MustRemove(tc.OutDir)
		}
	}
}
