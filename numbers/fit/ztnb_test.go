package fit

import (
	"math"
	"testing"
)

var ZeroTruncatedNegativeBinomialTests = []struct {
	R float64
	P float64
}{
	{
		R: 1.0,
		P: 0.4,
	},
	{
		R: 2,
		P: 0.1,
	},
	{
		R: 2,
		P: 0.4,
	},
	{
		R: 6,
		P: 0.4,
	},
	{
		R: 6,
		P: 0.4,
	},
}

func TestZeroTruncatedNegativeBinomial(t *testing.T) {
	var currR, currP float64
	var numVariates = 10000
	var countSlice []int
	var tmpCountSlice []int
	var currK int

	for _, v := range ZeroTruncatedNegativeBinomialTests {
		countSlice = make([]int, 1)
		countSlice[0] = 0
		for i := 0; i < numVariates; i++ {
			currK = randNegativeBinomial(v.R, v.P)
			if currK >= len(countSlice) {
				tmpCountSlice = make([]int, currK+1)
				copy(tmpCountSlice, countSlice)
				countSlice = tmpCountSlice
			}
			countSlice[currK]++
		}
		nbR, nbP, _ := NegativeBinomialFromCountSlice(countSlice)

		//var currLossPlot [][]float64
		//var err error
		_, _, lowestLossR, lowestLossP := plotLossSurfaceZTNB(countSlice, 0.1, 10, 0.1, 0.01, 0.99, 0.01)
		//out := fileio.EasyCreate(fmt.Sprintf(v.InFile + ".loss.txt"))

		//for i := range currLossPlot {
		//		for j := range currLossPlot[i] {
		//			_, err = fmt.Fprintf(out, "%v\t", currLossPlot[i][j])
		//			exception.PanicOnErr(err)
		//		}
		//		_, err = fmt.Fprintf(out, "\n")
		//	}

		//	err = out.Close()
		//	exception.PanicOnErr(err)

		currR, currP = ZeroTruncatedNegativeBinomial(countSlice, 5.0, 0.5, 0.1, 0.01)
		if math.Abs(currR-v.R) > 0.4 {
			t.Errorf("Error: NbR: %.1f. LowestLossR: %.1f. currR: %.1f. ExpectedR: %.1f.\n", nbR, lowestLossR, currR, v.R)
		}
		if math.Abs(currP-v.P) > 0.4 {
			t.Errorf("Error: NbP: %.1f. LowestLossP: %.1f. currP: %.1f. ExpectedP: %.1f.\n", nbP, lowestLossP, currP, v.P)
		}
	}
}
