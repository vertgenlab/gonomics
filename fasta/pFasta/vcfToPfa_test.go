package pFasta

import (
	"testing"
	"math/rand"
)

var vcfToPfaTests = []struct {
	InputVcf    string
	RefFa    string
	Start	int
	End	int
	SetSeed  int64
	Expected string
}{
	{InputVcf: "testdata/test_vcftoPfa_input_1.vcf",
		RefFa: "testdata/test_vcfToPfa_input_1.fa",
		Start: 1,
		End: 36,
		SetSeed: 7,
		Expected: "testdata/test_vcfToPfa_expected_1.vcf",
	},
	{InputVcf: "testdata/test_vcftoPfa_input_2.vcf",
		RefFa: "testdata/test_vcfToPfa_input_1.fa",
		Start: 1,
		End: 44,
		SetSeed: 7,
		Expected: "testdata/test_vcfToPfa_expected_2.vcf",
	},
}

func TestVcfToPfa(t *testing.T) {
	for _, tc := range vcfToPfaTests {
		rand.Seed(tc.SetSeed)
		observed := []PFasta{VcfToPfa(tc.InputVcf, tc.RefFa, tc.Start, tc.End)}
		Write("testdata/test_vcfToPfa_output_2.pfa", observed)
		if !fasta.IsEqual(observed, testCase.Expected) {
			t.Errorf("Error: in pFasta. Sample valid input test not as expected.\n")
		}
	}
}