package popgen

import (
	"github.com/vertgenlab/gonomics/numbers"
	"math"
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
		answer = numbers.AddLog(answer, fCache[j])
	}
	return answer
}

//AncestralAscertainmentDenominator is P(Asc | alpha), a constant normalizing factor for P(i | n, alpha) with the ancestral ascertainment correction.
func AncestralAscertainmentDenominator(fCache []float64, fCacheSum float64, d int) float64 {
	var answer float64 = math.Inf(-1)
	var current float64
	for j := 1; j < len(fCache); j++ {
		current = numbers.MultiplyLog(numbers.DivideLog(fCache[j], fCacheSum), AncestralAscertainmentProbability(len(fCache), j, d))
		answer = numbers.AddLog(answer, current)
	}
	return answer
}

//DerivedAscertainmentDenominator is P(Asc | alpha), a constant normalizing factor for P(i | n, alpha) with the derived ascertainment correction.
func DerivedAscertainmentDenominator(fCache []float64, fCacheSum float64, d int) float64 {
	var answer float64 = math.Inf(-1)
	var current float64
	for j := 1; j < len(fCache); j++ {
		current = numbers.MultiplyLog(numbers.DivideLog(fCache[j], fCacheSum), DerivedAscertainmentProbability(len(fCache), j, d))
		answer = numbers.AddLog(answer, current)
	}
	return answer
}

//AncestralAscertainmentProbability returns P(Asc | i, alpha) for ancestral allele ascertainment corrections.
func AncestralAscertainmentProbability(n int, i int, d int) float64 {
	return numbers.DivideLog(numbers.BinomCoefficientLog(n-i, d), numbers.BinomCoefficientLog(n, d))
}

//DerivedAscertainmentProbability returns P(Asc | i, alpha) for derived allele ascertainment corrections.
func DerivedAscertainmentProbability(n int, i int, d int) float64 {
	return numbers.DivideLog(numbers.BinomCoefficientLog(i, d), numbers.BinomCoefficientLog(n, d))
}

//AlleleFrequencyProbabilityAncestralAscertainment returns P(i | Asc, alpha) when the variant set has an ancestral allele ascertainment bias.
func AlleleFrequencyProbabilityAncestralAscertainment(alpha float64, i int, n int, d int, binomCache [][]float64, integralError float64) float64 {
	fCache := BuildFCache(n, alpha, binomCache, integralError)
	fCacheSum := GetFCacheSum(fCache)
	pIGivenAlpha := numbers.DivideLog(fCache[i], fCacheSum)

	return numbers.DivideLog(numbers.MultiplyLog(pIGivenAlpha, AncestralAscertainmentProbability(n, i, d)), AncestralAscertainmentDenominator(fCache, fCacheSum, d))
}

//AlleleFrequencyProbabilityDerivedAscertainment returns P(i | Asc, alpha) when the variant set has a derived allele ascertainment bias.
func AlleleFrequencyProbabilityDerivedAscertainment(alpha float64, i int, n int, d int, binomCache [][]float64, integralError float64) float64 {
	fCache := BuildFCache(n, alpha, binomCache, integralError)
	fCacheSum := GetFCacheSum(fCache)
	pIGivenAlpha := numbers.DivideLog(fCache[i], fCacheSum)

	return numbers.DivideLog(numbers.MultiplyLog(pIGivenAlpha, DerivedAscertainmentProbability(n, i, d)), DerivedAscertainmentDenominator(fCache, fCacheSum, d))
}

//AfsDivergenceAscertainmentLikelihood is like AfsLikelihood, but makes a correction for divergence-based ascertainment when variant sets were selected for divergence or identity between two groups of d individuals.
func AfsDivergenceAscertainmentLikelihood(afs Afs, alpha []float64, binomMap [][]float64, d int, integralError float64) float64 {
	var answer float64 = 0.0
	for j := 0; j < len(afs.Sites); j++ {
		var currLikelihood float64
		if afs.Sites[j].Divergent {
			currLikelihood = AlleleFrequencyProbabilityDerivedAscertainment(alpha[j], afs.Sites[j].I, afs.Sites[j].N, d, binomMap, integralError)
		} else {
			currLikelihood = AlleleFrequencyProbabilityAncestralAscertainment(alpha[j], afs.Sites[j].I, afs.Sites[j].N, d, binomMap, integralError)
		}
		answer = numbers.MultiplyLog(answer, currLikelihood)
	}
	return answer
}

//AfsDivergenceAscertainmentFixedAlpha returns the likelihood of observing a particular derived allele frequency spectrum for a given selection parameter alpha.
//This is the special case where every segregating site has the same value for selection. Also applies a correction for divergence-based ascertainment bias.
func AfsDivergenceAscertainmentFixedAlpha(afs Afs, alpha float64, binomMap [][]float64, d int, integralError float64) float64 {
	allN := findAllN(afs)
	var answer float64 = 0.0
	ancestralLikelihoodCache := BuildLikelihoodCache(allN)
	derivedLikelihoodCache := BuildLikelihoodCache(allN)
	for j := range afs.Sites {
		if afs.Sites[j].Divergent {
			if ancestralLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] == 0.0 {
				ancestralLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] = AlleleFrequencyProbabilityAncestralAscertainment(alpha, afs.Sites[j].I, afs.Sites[j].N, d, binomMap, integralError)
			}
			answer = numbers.MultiplyLog(answer, ancestralLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I])
		} else {
			if derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] == 0.0 {
				derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] = AlleleFrequencyProbabilityDerivedAscertainment(alpha, afs.Sites[j].I, afs.Sites[j].N, d, binomMap, integralError)
			}
			answer = numbers.MultiplyLog(answer, derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I])
		}
	}
	return answer
}

func AfsDerivedDivergenceAscertainmentFixedAlpha(afs Afs, alpha float64, binomMap [][]float64, d int, integralError float64) float64 {
	allN := findAllN(afs)
	var answer float64 = 0.0
	derivedLikelihoodCache := BuildLikelihoodCache(allN)
	uncorrectedCache := BuildLikelihoodCache(allN)
	for j := range afs.Sites {
		if afs.Sites[j].Divergent {
			if derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] == 0.0 {
				derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] = AlleleFrequencyProbabilityDerivedAscertainment(alpha, afs.Sites[j].I, afs.Sites[j].N, d, binomMap, integralError)
			}
			answer = numbers.MultiplyLog(answer, derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I])
		} else {
			if uncorrectedCache[afs.Sites[j].N][afs.Sites[j].I] == 0.0 {
				uncorrectedCache[afs.Sites[j].N][afs.Sites[j].I] = AlleleFrequencyProbability(afs.Sites[j].I, afs.Sites[j].N, alpha, binomMap, integralError)
			}
			answer = numbers.MultiplyLog(answer, uncorrectedCache[afs.Sites[j].N][afs.Sites[j].I])
		}
	}
	return answer
}

func AfsAncestralDivergenceAscertainmentFixedAlpha(afs Afs, alpha float64, binomMap [][]float64, d int, integralError float64) float64 {
	allN := findAllN(afs)
	var answer float64 = 0.0
	derivedLikelihoodCache := BuildLikelihoodCache(allN)
	uncorrectedCache := BuildLikelihoodCache(allN)
	for j := range afs.Sites {
		if !afs.Sites[j].Divergent { //site is not divergent, i.e. ref is ancestral state
			if derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] == 0.0 {
				derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I] = AlleleFrequencyProbabilityAncestralAscertainment(alpha, afs.Sites[j].I, afs.Sites[j].N, d, binomMap, integralError)
			}
			answer = numbers.MultiplyLog(answer, derivedLikelihoodCache[afs.Sites[j].N][afs.Sites[j].I])
		} else {
			if uncorrectedCache[afs.Sites[j].N][afs.Sites[j].I] == 0.0 {
				uncorrectedCache[afs.Sites[j].N][afs.Sites[j].I] = AlleleFrequencyProbability(afs.Sites[j].I, afs.Sites[j].N, alpha, binomMap, integralError)
			}
			answer = numbers.MultiplyLog(answer, uncorrectedCache[afs.Sites[j].N][afs.Sites[j].I])
		}
	}
	return answer
}

func AfsDerivedDivergenceAscertainmentFixedAlphaClosure(afs Afs, binomMap [][]float64, d int, integralError float64) func(float64) float64 {
	return func(alpha float64) float64 {
		return AfsDerivedDivergenceAscertainmentFixedAlpha(afs, alpha, binomMap, d, integralError)
	}
}

func AfsAncestralDivergenceAscertainmentFixedAlphaClosure(afs Afs, binomMap [][]float64, d int, integralError float64) func(float64) float64 {
	return func(alpha float64) float64 {
		return AfsAncestralDivergenceAscertainmentFixedAlpha(afs, alpha, binomMap, d, integralError)
	}
}

//AfsDivergenceAscertainmentFixedAlphaClosure returns a func(float64) float64 representing the likelihood function for a specific derived allele frequency spectrum with a single selection parameter alpha.
//Incorporates a site-by-site correction for divergence-based ascertainment bias.
func AfsDivergenceAscertainmentFixedAlphaClosure(afs Afs, binomMap [][]float64, d int, integralError float64) func(float64) float64 {
	return func(alpha float64) float64 {
		return AfsDivergenceAscertainmentFixedAlpha(afs, alpha, binomMap, d, integralError)
	}
}
