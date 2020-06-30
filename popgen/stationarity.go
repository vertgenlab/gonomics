package popgen

import (
	"math"
	"github.com/vertgenlab/gonomics/numbers"
)

/*
This file implements functions to construct a Hierarchical Bayesian model for inference of the selection parameter
alpha based on an allele frequency spectrum obtained from multiple alignment data.
More details can be obtained from the Ph.D. Thesis of Sol Katzman, available at https://compbio.soe.ucsc.edu/theses/Sol-Katzman-PhD-2010.pdf.
Equations from this thesis that are replicated here are marked.
*/

//eq 2.1
func AFSStationarity(p float64, alpha float64) float64 {
	return (1-math.Exp(-alpha * (1-p))) * 2 / ((1 - math.Exp(-alpha)) * p * (1 - p))
}

func AFSStationarityClosure(alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		return AFSStationarity(p, alpha)
	}
}

func AFSSampleClosure(n int, k int, alpha float64) func(float64) float64 {
	return func (p float64) float64 {
		return AFSStationarity(p, alpha) * numbers.BinomialDist(n, k, p)
	}
}

//eq. 2.2
func AFSSampleDensity(n int, k int, alpha float64) float64 {
	f := AFSSampleClosure(n, k, alpha)
	return numbers.DefiniteIntegral(f, 0, 1)
}

//eq 2.3
func AlleleFrequencyProbability(i int, n int, alpha float64) float64 {
	var denominator float64
	for j := 1; j < n - 1; j++ {
		denominator = denominator + AFSSampleDensity(j, n, alpha)
	}
	return AFSSampleDensity(i, n, alpha) / denominator
}


//eq 2.4
//afs array has a dummy variable in position 0, so loop starts at 1.
func AFSLikelihood(afs []int, alpha float64) float64 {
	var answer float64 = 1.0
	for i := 1; i < len(afs); i++ {
		answer = answer * AlleleFrequencyProbability(afs[i], len(afs), alpha)
	}
	return answer
}