package numbers

import (
	"fmt"
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
		variates = make([]int, v.VariateCount)
		for i = range v.ExpectedProb {
			if fmt.Sprintf("%.6g", currAlias.Probability[i]) != fmt.Sprintf("%.6g", v.ExpectedProb[i]) {
				t.Errorf("Error: RandBinomial produced an incorrect alias probability vector.\nExpected: %v.\nFound: %v.\n", v.ExpectedProb, currAlias.Probability)
			}
		}
		for i = range v.ExpectedAlias {
			if v.ExpectedAlias[i] != currAlias.Alias[i] {
				t.Errorf("Error: RandBinomial produced an incorrect alias index.\nExpected: %v.\nFound: %v.\n", v.ExpectedAlias, currAlias.Alias)
			}
		}

		// print variates for plotting
		for i = range variates {
			variates[i] = RandBinomial(currAlias)
			// fmt.Println(variates[i])
		}
	}
}
