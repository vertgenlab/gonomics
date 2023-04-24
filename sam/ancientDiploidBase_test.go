package sam

import (
	"math"
	"testing"
)

var AncientLikelihoodExpressionTests = []struct {
	Count                                                    int
	Epsilon                                                  float64
	Lambda                                                   float64
	ExpectedEpsilonOverThree                                 float64
	ExpectedOneMinusEpsilon                                  float64
	ExpectedOneMinusEpsilonMinusLambda                       float64
	ExpectedEpsilonOverThreePlusLambda                       float64
	ExpectedPointFiveMinusEpsilonOverThree                   float64
	ExpectedEpsilonOverThreePlusLambdaOverTwo                float64
	ExpectedPointFiveMinusEpsilonOverThreePlusLambdaOverTwo  float64
	ExpectedPointFiveMinusEpsilonOverThreeMinusLambdaOverTwo float64
}{
	{Count: 10,
		Epsilon:                                   0.01,
		Lambda:                                    0.05,
		ExpectedEpsilonOverThree:                  -57.03782474656201,
		ExpectedOneMinusEpsilon:                   -0.1005033585350145,
		ExpectedOneMinusEpsilonMinusLambda:        -0.6187540371808753,
		ExpectedEpsilonOverThreePlusLambda:        -29.311937524164197,
		ExpectedPointFiveMinusEpsilonOverThree:    -6.998361687107419,
		ExpectedEpsilonOverThreePlusLambdaOverTwo: -35.6371631115993,
		ExpectedPointFiveMinusEpsilonOverThreePlusLambdaOverTwo:  -6.507264646759933,
		ExpectedPointFiveMinusEpsilonOverThreeMinusLambdaOverTwo: -7.514827575729088,
	},
}

func TestAncientLikelihoodExpressions(t *testing.T) {
	var actual float64
	var cache = make([]float64, 0) // this will serve as an empty cache for all expression tests
	var ancientCache = AncientLikelihoodCache{
		EpsilonOverThree:                                 cache,
		OneMinusEpsilon:                                  cache,
		OneMinusEpsilonMinusLambda:                       cache,
		EpsilonOverThreePlusLambda:                       cache,
		PointFiveMinusEpsilonOverThree:                   cache,
		EpsilonOverThreePlusLambdaOverTwo:                cache,
		PointFiveMinusEpsilonOverThreePlusLambdaOverTwo:  cache,
		PointFiveMinusEpsilonOverThreeMinusLambdaOverTwo: cache,
	}
	for _, v := range AncientLikelihoodExpressionTests {
		actual = epsilonOverThreeLikelihoodExpression(v.Count, v.Epsilon, ancientCache)
		if math.Abs(actual-v.ExpectedEpsilonOverThree)/v.ExpectedEpsilonOverThree > 10e-6 {
			t.Errorf("Error: epsilonOverThreeLikelihoodExpression. Count: %v. Epsilon: %v. Found: %v. Expected: %v.\n", v.Count, v.Epsilon, actual, v.ExpectedEpsilonOverThree)
		}
		actual = oneMinusEpsilonLikelihoodExpression(v.Count, v.Epsilon, ancientCache)
		if math.Abs(actual-v.ExpectedOneMinusEpsilon)/v.ExpectedOneMinusEpsilon > 10e-6 {
			t.Errorf("Error: oneMinusEpsilonLikelihoodExpression. Count: %v. Epsilon: %v. Found: %v. Expected: %v.\n", v.Count, v.Epsilon, actual, v.ExpectedOneMinusEpsilon)
		}
		actual = oneMinusEpsilonMinusLambdaLikelihoodExpression(v.Count, v.Epsilon, v.Lambda, ancientCache)
		if math.Abs(actual-v.ExpectedOneMinusEpsilonMinusLambda)/v.ExpectedOneMinusEpsilonMinusLambda > 10e-6 {
			t.Errorf("Error: oneMinusEpsilonMinusLambdaLikelihoodExpression. Count: %v. Epsilon: %v. Lambda: %v. Found: %v. Expected: %v.\n", v.Count, v.Epsilon, v.Lambda, actual, v.ExpectedOneMinusEpsilonMinusLambda)
		}
		actual = epsilonOverThreePlusLambdaLikelihoodExpression(v.Count, v.Epsilon, v.Lambda, ancientCache)
		if math.Abs(actual-v.ExpectedEpsilonOverThreePlusLambda)/v.ExpectedEpsilonOverThreePlusLambda > 10e-6 {
			t.Errorf("Error: epsilonOverThreePlusLambdaLikelihoodExpression. Count: %v. Epsilon: %v. Lambda: %v. Found: %v. Expected: %v.\n", v.Count, v.Epsilon, v.Lambda, actual, v.ExpectedEpsilonOverThreePlusLambda)
		}
		actual = pointFiveMinusEpsilonOverThreeLikelihoodExpression(v.Count, v.Epsilon, ancientCache)
		if math.Abs(actual-v.ExpectedPointFiveMinusEpsilonOverThree)/v.ExpectedPointFiveMinusEpsilonOverThree > 10e-6 {
			t.Errorf("Error: pointFiveMinusEpsilonOverThreeLikelihoodExpression. Count: %v. Epsilon: %v. Found: %v. Expected: %v.\n", v.Count, v.Epsilon, actual, v.ExpectedPointFiveMinusEpsilonOverThree)
		}
		actual = epsilonOverThreePlusLambdaOverTwoLikelihoodExpression(v.Count, v.Epsilon, v.Lambda, ancientCache)
		if math.Abs(actual-v.ExpectedEpsilonOverThreePlusLambdaOverTwo)/v.ExpectedEpsilonOverThreePlusLambdaOverTwo > 10e-6 {
			t.Errorf("Error: epsilonOverThreePlusLambdaOverTwoLikelihoodExpression. Count: %v. Epsilon: %v. Lambda: %v. Found: %v. Expected: %v.\n", v.Count, v.Epsilon, v.Lambda, actual, v.ExpectedEpsilonOverThreePlusLambdaOverTwo)
		}
		actual = pointFiveMinusEpsilonOverThreePlusLambdaOverTwoLikelihoodExpression(v.Count, v.Epsilon, v.Lambda, ancientCache)
		if math.Abs(actual-v.ExpectedPointFiveMinusEpsilonOverThreePlusLambdaOverTwo)/v.ExpectedPointFiveMinusEpsilonOverThreePlusLambdaOverTwo > 10e-6 {
			t.Errorf("Error: pointFiveMinusEpsilonOverThreePlusLambdaOverTwoLikelihoodExpression. Count: %v. Epsilon: %v. Lambda: %v. Found: %v. Expected: %v.\n", v.Count, v.Epsilon, v.Lambda, actual, v.ExpectedPointFiveMinusEpsilonOverThreePlusLambdaOverTwo)
		}
		actual = pointFiveMinusEpsilonOverThreeMinusLambdaOverTwoLikelihoodExpression(v.Count, v.Epsilon, v.Lambda, ancientCache)
		if math.Abs(actual-v.ExpectedPointFiveMinusEpsilonOverThreeMinusLambdaOverTwo)/v.ExpectedPointFiveMinusEpsilonOverThreeMinusLambdaOverTwo > 10e-6 {
			t.Errorf("Error: pointFiveMinusEpsilonOverThreeMinusLambdaOverTwoLikelihoodExpression. Count: %v. Epsilon: %v. Lambda: %v. Found: %v. Expected: %v.\n", v.Count, v.Epsilon, v.Lambda, actual, v.ExpectedPointFiveMinusEpsilonOverThreeMinusLambdaOverTwo)
		}
	}
}
