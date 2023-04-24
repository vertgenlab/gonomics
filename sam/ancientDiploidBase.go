package sam

import (
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"log"
	"math"
)

// AncientLikelihoodCache stores the following short expressions from the likelihood functions to
// save on repeated logspace.Pow calls.
type AncientLikelihoodCache struct {
	EpsilonOverThree                                 []float64 // \frac{\epsilon}{3}
	OneMinusEpsilon                                  []float64 // 1 - \epsilon
	OneMinusEpsilonMinusLambda                       []float64 // 1 - \frac{\epsilon}{3} - \lambda
	EpsilonOverThreePlusLambda                       []float64 // \frac{\epsilon}{3} + \lambda
	PointFiveMinusEpsilonOverThree                   []float64 // 0.5 - \frac{\epsilon}{3}
	EpsilonOverThreePlusLambdaOverTwo                []float64 // \frac{\epsilon}{3} + \frac{\lambda}{2}
	PointFiveMinusEpsilonOverThreePlusLambdaOverTwo  []float64 // 0.5 - \frac{\epsilon}{3} + \frac{\lambda}{2}
	PointFiveMinusEpsilonOverThreeMinusLambdaOverTwo []float64 // 0.5 - \frac{\epsilon}{3} - \frac{\lambda}{2}
}

// ancientBaseLikelihood is a helper function of DiploidBaseCallFromPile. For a given genotype (geno), and counts
// for the four bases, a PCR/sequencing error rate epsilon, and a cytosine deamination rate lambda, this function
// calculates the genotype likelihood.
// Individual terms of the likelihood expression are calculated using an AncientLikelihoodCache struct to save on
// repeated logspace.Pow calls.
func ancientBaseLikelihood(aCount int, cCount int, gCount int, tCount int, geno DiploidBase, epsilon float64, lambda float64, cache AncientLikelihoodCache) float64 {
	var firstTerm, secondTerm, thirdTerm, fourthTerm float64 = 0, 0, 0, 0 //default to 1 in logspace
	switch geno {
	case AA:
		firstTerm = epsilonOverThreeLikelihoodExpression(cCount+gCount+tCount, epsilon, cache)
		secondTerm = oneMinusEpsilonLikelihoodExpression(aCount, epsilon, cache)
	case AC:
		firstTerm = pointFiveMinusEpsilonOverThreeLikelihoodExpression(aCount, epsilon, cache)
		secondTerm = pointFiveMinusEpsilonOverThreeMinusLambdaOverTwoLikelihoodExpression(cCount, epsilon, lambda, cache)
		thirdTerm = epsilonOverThreeLikelihoodExpression(gCount, epsilon, cache)
		fourthTerm = epsilonOverThreePlusLambdaOverTwoLikelihoodExpression(tCount, epsilon, lambda, cache)
	case AG:
		firstTerm = pointFiveMinusEpsilonOverThreePlusLambdaOverTwoLikelihoodExpression(aCount, epsilon, lambda, cache)
		secondTerm = epsilonOverThreeLikelihoodExpression(cCount+gCount, epsilon, cache)
		thirdTerm = pointFiveMinusEpsilonOverThreeMinusLambdaOverTwoLikelihoodExpression(gCount, epsilon, lambda, cache)
	case AT:
		firstTerm = pointFiveMinusEpsilonOverThreeLikelihoodExpression(aCount+tCount, epsilon, cache)
		secondTerm = epsilonOverThreeLikelihoodExpression(cCount+gCount, epsilon, cache)
	case CC:
		firstTerm = epsilonOverThreeLikelihoodExpression(aCount+gCount, epsilon, cache)
		secondTerm = oneMinusEpsilonMinusLambdaLikelihoodExpression(cCount, epsilon, lambda, cache)
		thirdTerm = epsilonOverThreePlusLambdaLikelihoodExpression(tCount, epsilon, lambda, cache)
	case CG:
		firstTerm = epsilonOverThreePlusLambdaOverTwoLikelihoodExpression(aCount, epsilon, lambda, cache)
		secondTerm = pointFiveMinusEpsilonOverThreeMinusLambdaOverTwoLikelihoodExpression(cCount+gCount, epsilon, lambda, cache)
		thirdTerm = epsilonOverThreePlusLambdaOverTwoLikelihoodExpression(tCount, epsilon, lambda, cache)
	case CT:
		firstTerm = epsilonOverThreeLikelihoodExpression(aCount+gCount, epsilon, cache)
		secondTerm = pointFiveMinusEpsilonOverThreeMinusLambdaOverTwoLikelihoodExpression(cCount, epsilon, lambda, cache)
		thirdTerm = pointFiveMinusEpsilonOverThreePlusLambdaOverTwoLikelihoodExpression(tCount, epsilon, lambda, cache)
	case GG:
		firstTerm = epsilonOverThreePlusLambdaLikelihoodExpression(aCount, epsilon, lambda, cache)
		secondTerm = epsilonOverThreeLikelihoodExpression(cCount+tCount, epsilon, cache)
		thirdTerm = oneMinusEpsilonMinusLambdaLikelihoodExpression(gCount, epsilon, lambda, cache)
	case GT:
		firstTerm = epsilonOverThreePlusLambdaOverTwoLikelihoodExpression(aCount, epsilon, lambda, cache)
		secondTerm = epsilonOverThreeLikelihoodExpression(cCount, epsilon, cache)
		thirdTerm = pointFiveMinusEpsilonOverThreeMinusLambdaOverTwoLikelihoodExpression(gCount, epsilon, lambda, cache)
		fourthTerm = pointFiveMinusEpsilonOverThreeLikelihoodExpression(tCount, epsilon, cache)
	case TT:
		firstTerm = epsilonOverThreeLikelihoodExpression(aCount+cCount+gCount, epsilon, cache)
		secondTerm = oneMinusEpsilonLikelihoodExpression(tCount, epsilon, cache)
	default:
		log.Fatalf("Error: Unknown genotype: %v.\n", geno)
	}
	return logspace.Multiply(firstTerm, logspace.Multiply(secondTerm, logspace.Multiply(thirdTerm, fourthTerm)))
}

// epsilonOverThreeLikelihoodExpression is a helper function of ancientBaseLikelihood that calculates the expression
// (epsilon/3.0)^count.
// Checks the cache to see if this expression has already been calculated. Otherwise, updates the cache.
func epsilonOverThreeLikelihoodExpression(count int, epsilon float64, cache AncientLikelihoodCache) float64 {
	if count < len(cache.EpsilonOverThree) { // if the coverage is within the cache bounds
		if cache.EpsilonOverThree[count] != 0 {
			return cache.EpsilonOverThree[count]
		} else {
			cache.EpsilonOverThree[count] = logspace.Pow(math.Log(epsilon/3.0), float64(count))
			return cache.EpsilonOverThree[count]
		}
	} else { // otherwise, calculate by hand
		return logspace.Pow(math.Log(epsilon/3.0), float64(count))
	}
}

// oneMinusEpsilonLikelihoodExpression is a helper function of ancientBaseLikelihood that calculates the expression
// (1-epsilon)^count.
// Checks the cache to see if this expression has already been calculated. Otherwise, updates the cache.
func oneMinusEpsilonLikelihoodExpression(count int, epsilon float64, cache AncientLikelihoodCache) float64 {
	if count < len(cache.OneMinusEpsilon) { // if the coverage is within the cache bounds
		if cache.OneMinusEpsilon[count] != 0 {
			return cache.OneMinusEpsilon[count]
		} else {
			cache.OneMinusEpsilon[count] = logspace.Pow(math.Log(1.0-epsilon), float64(count))
			return cache.OneMinusEpsilon[count]
		}
	} else {
		return logspace.Pow(math.Log(1.0-epsilon), float64(count))
	}
}

// oneMinusEpsilonMinusLambdaLikelihoodExpression is a helper function of ancientBaseLikelihood that calculates the expression
// (1-epsilon-lambda)^count.
// Checks the cache to see if this expression has already been calculated. Otherwise, updates the cache.
func oneMinusEpsilonMinusLambdaLikelihoodExpression(count int, epsilon float64, lambda float64, cache AncientLikelihoodCache) float64 {
	if count < len(cache.OneMinusEpsilonMinusLambda) { // if the coverage is within the cache bounds
		if cache.OneMinusEpsilonMinusLambda[count] != 0 {
			return cache.OneMinusEpsilonMinusLambda[count]
		} else {
			cache.OneMinusEpsilonMinusLambda[count] = logspace.Pow(math.Log(1.0-epsilon-lambda), float64(count))
			return cache.OneMinusEpsilonMinusLambda[count]
		}
	} else {
		return logspace.Pow(math.Log(1.0-epsilon-lambda), float64(count))
	}
}

// epsilonOverThreePlusLambdaLikelihoodExpression is a helper function of ancientBaseLikelihood that calculates the expression
// (epsilon/3 + lambda)^count.
// Checks the cache to see if this expression has already been calculated. Otherwise, updates the cache.
func epsilonOverThreePlusLambdaLikelihoodExpression(count int, epsilon float64, lambda float64, cache AncientLikelihoodCache) float64 {
	if count < len(cache.EpsilonOverThreePlusLambda) { // if the coverage is within the cache bounds
		if cache.EpsilonOverThreePlusLambda[count] != 0 {
			return cache.EpsilonOverThreePlusLambda[count]
		} else {
			cache.EpsilonOverThreePlusLambda[count] = logspace.Pow(math.Log((epsilon/3.0)+lambda), float64(count))
			return cache.EpsilonOverThreePlusLambda[count]
		}
	} else {
		return logspace.Pow(math.Log((epsilon/3.0)+lambda), float64(count))
	}
}

// pointFiveMinusEpsilonOverThreeLikelihoodExpression is a helper function of ancientBaseLikelihood that calculates the expression
// (0.5 - epsilon/3.0)^count.
// Checks the cache to see if this expression has already been calculated. Otherwise, updates the cache.
func pointFiveMinusEpsilonOverThreeLikelihoodExpression(count int, epsilon float64, cache AncientLikelihoodCache) float64 {
	if count < len(cache.PointFiveMinusEpsilonOverThree) { // if the coverage is within the cache bounds
		if cache.PointFiveMinusEpsilonOverThree[count] != 0 {
			return cache.PointFiveMinusEpsilonOverThree[count]
		} else {
			cache.PointFiveMinusEpsilonOverThree[count] = logspace.Pow(math.Log(0.5-(epsilon/3.0)), float64(count))
			return cache.PointFiveMinusEpsilonOverThree[count]
		}
	} else {
		return logspace.Pow(math.Log(0.5-(epsilon/3.0)), float64(count))
	}
}

// epsilonOverThreePlusLambdaOverTwoLikelihoodExpression is a helper function of ancientBaseLikelihood that calculates the expression
// (epsilon/3.0 + lambda/2.0)^count.
// Checks the cache to see if this expression has already been calculated. Otherwise, updates the cache.
func epsilonOverThreePlusLambdaOverTwoLikelihoodExpression(count int, epsilon float64, lambda float64, cache AncientLikelihoodCache) float64 {
	if count < len(cache.EpsilonOverThreePlusLambdaOverTwo) {
		if cache.EpsilonOverThreePlusLambdaOverTwo[count] != 0 {
			return cache.EpsilonOverThreePlusLambdaOverTwo[count]
		} else {
			cache.EpsilonOverThreePlusLambdaOverTwo[count] = logspace.Pow(math.Log((epsilon/3.0)+(lambda/2.0)), float64(count))
			return cache.EpsilonOverThreePlusLambdaOverTwo[count]
		}
	} else {
		return logspace.Pow(math.Log((epsilon/3.0)+(lambda/2.0)), float64(count))
	}
}

// pointFiveMinusEpsilonOverThreePlusLambdaOverTwoLikelihoodExpression is a helper function of ancientBaseLikelihood that calculates the expression
// (0.5 - epsilon/3.0 + lambda/2.0)^count.
// Checks the cache to see if this expression has already been calculated. Otherwise, updates the cache.
func pointFiveMinusEpsilonOverThreePlusLambdaOverTwoLikelihoodExpression(count int, epsilon float64, lambda float64, cache AncientLikelihoodCache) float64 {
	if count < len(cache.PointFiveMinusEpsilonOverThreePlusLambdaOverTwo) {
		if cache.PointFiveMinusEpsilonOverThreePlusLambdaOverTwo[count] != 0 {
			return cache.PointFiveMinusEpsilonOverThreePlusLambdaOverTwo[count]
		} else {
			cache.PointFiveMinusEpsilonOverThreePlusLambdaOverTwo[count] = logspace.Pow(math.Log(0.5-(epsilon/3.0)+(lambda/2.0)), float64(count))
			return cache.PointFiveMinusEpsilonOverThreePlusLambdaOverTwo[count]
		}
	} else {
		return logspace.Pow(math.Log(0.5-(epsilon/3.0)+(lambda/2.0)), float64(count))
	}
}

// pointFiveMinusEpsilonOverThreeMinusLambdaOverTwoLikelihoodExpression is a helper function of ancientBaseLikelihood that calculates the expression
// (0.5 - epsilon/3.0 - lambda/2.0)^count.
// Checks the cache to see if this expression has already been calculated. Otherwise, updates the cache.
func pointFiveMinusEpsilonOverThreeMinusLambdaOverTwoLikelihoodExpression(count int, epsilon float64, lambda float64, cache AncientLikelihoodCache) float64 {
	if count < len(cache.PointFiveMinusEpsilonOverThreeMinusLambdaOverTwo) {
		if cache.PointFiveMinusEpsilonOverThreeMinusLambdaOverTwo[count] != 0 {
			return cache.PointFiveMinusEpsilonOverThreeMinusLambdaOverTwo[count]
		} else {
			cache.PointFiveMinusEpsilonOverThreeMinusLambdaOverTwo[count] = logspace.Pow(math.Log(0.5-(epsilon/3.0)-(lambda/2.0)), float64(count))
			return cache.PointFiveMinusEpsilonOverThreeMinusLambdaOverTwo[count]
		}
	} else {
		return logspace.Pow(math.Log(0.5-(epsilon/3.0)-(lambda/2.0)), float64(count))
	}
}
