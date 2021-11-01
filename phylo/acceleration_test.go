package phylo

import "testing"

var BranchLengthsAlternatingLeastSquaresTests = []struct {
	D                     AccelDistancesAndWeights
	AllowNegative         bool
	Epsilon               float64
	ExpectedBranchLengths AccelBranchLengths
}{
	{AccelDistancesAndWeights{ //This particular test was one that on genome-wide data failed to converge at 1000 iterations. Slowly descends on minimum, achieves in 10k iterations.
		62,
		67,
		37,
		31,
		63,
		71,
		calculateWeight(62, 1000),
		calculateWeight(67, 1000),
		calculateWeight(37, 1000),
		calculateWeight(31, 1000),
		calculateWeight(63, 1000),
		calculateWeight(71, 1000),
	}, false, 1e-8, AccelBranchLengths{24.179037634841176, 19.55544803664856, 0, 18.710259085171323, 21.978125070026614}},
}

func TestBranchLengthsAlternatingLeastSquares(t *testing.T) {
	var currAnswer AccelBranchLengths
	for _, v := range BranchLengthsAlternatingLeastSquaresTests {
		currAnswer = BranchLengthsAlternatingLeastSquares(v.D, v.AllowNegative, false, 1000, v.Epsilon)
		if currAnswer != v.ExpectedBranchLengths {
			t.Errorf("Error in BranchLengthsAlternatingLeastSquares.\nAnswer: %v.\nExpected: %v.", currAnswer, v.ExpectedBranchLengths)
		}
	}
}
