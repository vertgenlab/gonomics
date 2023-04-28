package sam

import (
	"math"
	"testing"
)

var AncientBaseLikelihoodTests = []struct {
	aCount         int
	cCount         int
	gCount         int
	tCount         int
	epsilon        float64
	lambda         float64
	ExpectedValues []float64
}{
	{
		aCount:         13,
		cCount:         5,
		gCount:         0,
		tCount:         2,
		epsilon:        0.01,
		lambda:         0.05,
		ExpectedValues: []float64{-40.057131688688926, -19.982716603424052, -48.38592136338132, -39.01645490394213, -80.3209366939539, -57.213158455263496, -79.20803888774715, -78.03199610400686, -76.24689675578159, -102.68818521551862},
	},
}

func TestAncientBaseLikelihood(t *testing.T) {
	var actual float64
	var cache = make([]float64, 0)
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
	for _, v := range AncientBaseLikelihoodTests {
		for geno := 0; geno < 10; geno++ {
			actual = ancientBaseLikelihood(v.aCount, v.cCount, v.gCount, v.tCount, DiploidBase(geno), v.epsilon, v.lambda, ancientCache)
			if math.Abs(actual-v.ExpectedValues[geno])/v.ExpectedValues[geno] > 1e-6 {
				t.Errorf("Error: ancientBaseLikelihood. Expected: %v. Found: %v.\n", v.ExpectedValues[geno], actual)
			}
		}
	}
}

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
