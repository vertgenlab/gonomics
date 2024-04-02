package simulate

import (
	"github.com/vertgenlab/gonomics/numbers/matrix"
	"testing"
)

var ParseSubstitutionMatrixTests = []struct {
	InFile    string
	Expected  [][]float64
	Precision float64
}{
	{InFile: "testdata/substitutionMatrix.txt",
		Expected: [][]float64{
			{0.90, 0.05, 0.02, 0.03},
			{0.05, 0.90, 0.03, 0.02},
			{0.02, 0.03, 0.90, 0.05},
			{0.03, 0.02, 0.05, 0.90},
		},
		Precision: 1e-3,
	},
}

func TestParseSubstitutionMatrix(t *testing.T) {
	var observed [][]float64
	for _, v := range ParseSubstitutionMatrixTests {
		observed = ParseSubstitutionMatrix(v.InFile)
		if !matrix.ApproxEqual(observed, v.Expected, v.Precision) {
			t.Errorf("Error: parseSubstitutionMatrix. Output was not as expected.\n%v\n", observed)
		}
	}
}
