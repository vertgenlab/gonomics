package numbers

import (
	"testing"
)

var RandBinomialTests = []struct {
	N             int
	P             float64
	VariateCount  int
	ExpectedProb  []float64
	ExpectedAlias []int
}{
	{N: 5,
		P:             0.1,
		VariateCount:  20,
		ExpectedProb:  []float64{1, 0.9683600000000008, 0.4374000000000003, 0.04860000000000002, 0.0027000000000000014, 6.0000000000000096e-05},
		ExpectedAlias: []int{0, 0, 0, 0, 0, 1},
	},
}

func TestRandBinomial(t *testing.T) {
	var currAlias BinomialAlias
	var variates []int
	var i int
	for _, v := range RandBinomialTests {
		currAlias = MakeBinomialAlias(v.N, v.P)
		// Generate a large number of samples
		variates = make([]int, v.VariateCount)
		// Check if the number of elements in Probability and Alias slices are correct
		if len(currAlias.Probability) != v.N+1 {
			t.Errorf("Error: RandBinomial produced incorrect number of probabilities in the alias probability vector.\nExpected: %d.\nFound: %d.\n", v.N+1, len(currAlias.Probability))
		}
		if len(currAlias.Alias) != v.N+1 {
			t.Errorf("Error: RandBinomial produced incorrect number of aliases.\nExpected: %d.\nFound: %d.\n", v.N+1, len(currAlias.Alias))
		}
		// Check probabilities and aliases
		for i = 0; i < v.N+1; i++ {
			if !ApproxEqual(currAlias.Probability[i], v.ExpectedProb[i], 0.000001) { // Use a tolerance for floating-point comparison
				t.Errorf("Error: RandBinomial produced an incorrect alias probability at index %d.\nExpected: %v.\nFound: %v.\n", i, v.ExpectedProb[i], currAlias.Probability[i])
			}
			if v.ExpectedAlias[i] != currAlias.Alias[i] {
				t.Errorf("Error: RandBinomial produced an incorrect alias index at index %d.\nExpected: %v.\nFound: %v.\n", i, v.ExpectedAlias[i], currAlias.Alias[i])
			}
		}

		for i = 0; i < len(variates); i++ {
			variates[i] = RandBinomial(currAlias)
			if variates[i] < 0 || variates[i] > v.N {
				t.Errorf("Error: RandBinomial produced a variate outside the valid range [0, %d]: %d\n", v.N, variates[i])
			}

		}
	}
}
