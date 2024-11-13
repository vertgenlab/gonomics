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
	var s VcfToPfaSettings
	for _, tc := range vcfToPfaTests {
		s = VcfToPfaSettings{
			InFile: tc.InFile,
			RefFile: tc.RefFile,
			OutDir: tc.OutDir,
			Start: tc.Start,
			End: tc.End,
		}

		vcfToPfa(s)
		observed := pFasta.Read(tc.OutDir)
		expected := pFasta.Read(tc.Expected)
		if !pFasta.AllAreEqual(observed, expected, tc.Precision) {
			t.Errorf("Error: in pFasta. Sample valid input test not as expected.\n")
		} else {
			fileio.MustRemove(tc.OutDir)
		}
	}
}
