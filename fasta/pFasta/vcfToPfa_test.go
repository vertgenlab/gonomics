package pFasta

import (
	"testing"
	"math/rand"
	"github.com/vertgenlab/gonomics/fileio"
)

var vcfToPfaTests = []struct {
	InputVcf    string
	RefFa    string
	Start	int
	End	int
	SetSeed  int64
	Expected string
	Output string
	Precision	float32
}{
	{InputVcf: "testdata/test_vcftoPfa_input_1.vcf",
		RefFa: "testdata/test_vcfToPfa_input_1.fa",
		Start: 1,
		End: 36,
		SetSeed: 7,
		Expected: "testdata/test_vcfToPfa_expected_1.pfa",
		Output: "testdata/test_vcfToPfa_observed_1.pfa",
		Precision: 1e-3,
	},
	{InputVcf: "testdata/test_vcftoPfa_input_2.vcf",
		RefFa: "testdata/test_vcfToPfa_input_1.fa",
		Start: 1,
		End: 44,
		SetSeed: 7,
		Expected: "testdata/test_vcfToPfa_expected_2.pfa",
		Output: "testdata/test_vcfToPfa_observed_2.pfa",
		Precision: 1e-3,
	},
}

func TestVcfToPfa(t *testing.T) {
	for _, tc := range vcfToPfaTests {
		rand.Seed(tc.SetSeed)
		observed := []PFasta{VcfToPfa(tc.InputVcf, tc.RefFa, tc.Start, tc.End)}
		Write(tc.Output, observed)
		expected := Read(tc.Expected)
		if !AllAreEqual(observed, expected, tc.Precision) {
			t.Errorf("Error: in pFasta. Sample valid input test not as expected.\n")
		} else {
			fileio.MustRemove(tc.Output)
		}
	}
}