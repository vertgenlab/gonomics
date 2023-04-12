package popgen

import (
	"log"
	"math"

	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
)

//BuildFCache builds a slice of len(n) where each index i contains log(F(i | n, alpha)), where F is popgen.AFSSampleDensity
func BuildFCache(n int, alpha float64, binomCache [][]float64, integralError float64) []float64 {
	var answer []float64 = make([]float64, n, n)
	for j := 1; j < n; j++ {
		answer[j] = AfsSampleDensity(n, j, alpha, binomCache, integralError)
	}
	return answer
}

//GetFCacheSum uses the fCAche built in BuildFCache and calculates the sum from j=1 to n-1
func GetFCacheSum(fCache []float64) float64 {
	var answer float64 = math.Inf(-1)
	for j := 1; j < len(fCache); j++ {
		answer = logspace.Add(answer, fCache[j])
	}
	return answer
}

//AncestralAscertainmentDenominator is P(Asc | alpha), a constant normalizing factor for P(i | n, alpha) with the ancestral ascertainment correction.
func AncestralAscertainmentDenominator(fCache []float64, fCacheSum float64, d int) float64 {
	var answer float64 = math.Inf(-1)
	var current float64
	for j := 1; j < len(fCache); j++ {
		current = logspace.Multiply(logspace.Divide(fCache[j], fCacheSum), AncestralAscertainmentProbability(len(fCache), j, d))
		answer = logspace.Add(answer, current)
	}
	return answer
}

//DerivedAscertainmentDenominator is P(Asc | alpha), a constant normalizing factor for P(i | n, alpha) with the derived ascertainment correction.
func DerivedAscertainmentDenominator(fCache []float64, fCacheSum float64, d int) float64 {
	var answer float64 = math.Inf(-1)
	var current float64
	for j := 1; j < len(fCache); j++ {
		current = logspace.Multiply(logspace.Divide(fCache[j], fCacheSum), DerivedAscertainmentProbability(len(fCache), j, d))
		answer = logspace.Add(answer, current)
	}
	return answer
}

//AncestralAscertainmentProbability returns P(Asc | i, alpha) for ancestral allele ascertainment corrections.
func AncestralAscertainmentProbability(n int, i int, d int) float64 {
	return logspace.Divide(numbers.BinomCoefficientLog(n-i, d), numbers.BinomCoefficientLog(n, d))
}

//DerivedAscertainmentProbability returns P(Asc | i, alpha) for derived allele ascertainment corrections.
func DerivedAscertainmentProbability(n int, i int, d int) float64 {
	return logspace.Divide(numbers.BinomCoefficientLog(i, d), numbers.BinomCoefficientLog(n, d))
}

//AlleleFrequencyProbabilityAncestralAscertainment returns P(i | Asc, alpha) when the variant set has an ancestral allele ascertainment bias.
func AlleleFrequencyProbabilityAncestralAscertainment(alpha float64, i int, n int, d int, binomCache [][]float64, integralError float64) float64 {
	fCache := BuildFCache(n, alpha, binomCache, integralError)
	fCacheSum := GetFCacheSum(fCache)
	pIGivenAlpha := logspace.Divide(fCache[i], fCacheSum)

	return logspace.Divide(logspace.Multiply(pIGivenAlpha, AncestralAscertainmentProbability(n, i, d)), AncestralAscertainmentDenominator(fCache, fCacheSum, d))
}

//AlleleFrequencyProbabilityDerivedAscertainment returns P(i | Asc, alpha) when the variant set has a derived allele ascertainment bias.
func AlleleFrequencyProbabilityDerivedAscertainment(alpha float64, i int, n int, d int, binomCache [][]float64, integralError float64) float64 {
	fCache := BuildFCache(n, alpha, binomCache, integralError)
	fCacheSum := GetFCacheSum(fCache)
	pIGivenAlpha := logspace.Divide(fCache[i], fCacheSum)

	return logspace.Divide(logspace.Multiply(pIGivenAlpha, DerivedAscertainmentProbability(n, i, d)), DerivedAscertainmentDenominator(fCache, fCacheSum, d))
}

//AfsDivergenceAscertainmentLikelihood is like AfsLikelihood, but makes a correction for divergence-based ascertainment when variant sets were selected for divergence or identity between two groups of d individuals.
func AfsDivergenceAscertainmentLikelihood(afs Afs, alpha []float64, binomMap [][]float64, d int, integralError float64) float64 {
	var answer float64 = 0.0
	for j := 0; j < len(afs.Sites); j++ {
		var currLikelihood float64
		if afs.Sites[j].L == Uncorrected {
			currLikelihood = AlleleFrequencyProbability(afs.Sites[j].I, afs.Sites[j].N, alpha[j], binomMap, integralError)
		} else if afs.Sites[j].L == Ancestral {
			currLikelihood = AlleleFrequencyProbabilityAncestralAscertainment(alpha[j], afs.Sites[j].I, afs.Sites[j].N, d, binomMap, integralError)
		} else if afs.Sites[j].L == Derived {
			currLikelihood = AlleleFrequencyProbabilityDerivedAscertainment(alpha[j], afs.Sites[j].I, afs.Sites[j].N, d, binomMap, integralError)
		} else {
			log.Fatalf("invalid likelihood field in segregating site struct")
		}
		answer = logspace.Multiply(answer, currLikelihood)
	}
	return answer
}

//AfsDivergenceAscertainmentFixedAlpha returns the likelihood of observing a particular derived allele frequency spectrum for a given selection parameter alpha.
//This is the special case where every segregating site has the same value for selection. Also applies a correction for divergence-based ascertainment bias.
func AfsDivergenceAscertainmentFixedAlpha(afs Afs, alpha float64, binomMap [][]float64, d int, integralError float64) float64 {
	allN := findAllN(afs)
	var answer float64 = 0.0
	uncorrectedLikelihoodCache := BuildLikelihoodCache(allN)
	ancestralLikelihoodCache := BuildLikelihoodCache(allN)
	derivedLikelihoodCache := BuildLikelihoodCache(allN)
	for j := range afs.Sites {
		switch afs.Sites[j].L {
		case Uncorrected:
			if uncorrectedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] == 0.0 { //if this particular segregating site has not already had its likelihood value cached, we want to calculate and cache it.
				uncorrectedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] = AlleleFrequencyProbability(afs.Sites[j].I, afs.Sites[j].N, alpha, binomMap, integralError)
			}
			answer = logspace.Multiply(answer, uncorrectedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I])
		case Ancestral:
			if ancestralLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] == 0.0 {
				ancestralLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] = AlleleFrequencyProbabilityAncestralAscertainment(alpha, afs.Sites[j].I, afs.Sites[j].N, d, binomMap, integralError)
			}
			answer = logspace.Multiply(answer, ancestralLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I])
		case Derived:
			if derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] == 0.0 {
				derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] = AlleleFrequencyProbabilityDerivedAscertainment(alpha, afs.Sites[j].I, afs.Sites[j].N, d, binomMap, integralError)
			}
			answer = logspace.Multiply(answer, derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I])
		}
	}
	return answer
}

//AfsDivergenceAscertainmentFixedAlphaClosure returns a func(float64) float64 representing the likelihood function for a specific derived allele frequency spectrum with a single selection parameter alpha.
//Incorporates a site-by-site correction for divergence-based ascertainment bias.
func AfsDivergenceAscertainmentFixedAlphaClosure(afs Afs, binomMap [][]float64, d int, integralError float64) func(float64) float64 {
	return func(alpha float64) float64 {
		return AfsDivergenceAscertainmentFixedAlpha(afs, alpha, binomMap, d, integralError)
	}
}
