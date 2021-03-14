package popgen

import "github.com/vertgenlab/gonomics/numbers"

func AncestralAscertainmentProbability(n int, i int, d int) float64 {
	return numbers.DivideLog(numbers.BinomCoefficientLog(n-i, d), numbers.BinomCoefficientLog(n, d))
}

func DerivedAscertainmentProbability(n int, i int, d int) float64 {
	return numbers.DivideLog(numbers.BinomCoefficientLog(i, d), numbers.BinomCoefficientLog(n, d))
}

func AncestralAscertainmentDenominator(alpha float64, n int, d int, binomCache [][]float64) float64 {
	var answer float64 = 1.0 //0 in logSpace
	var iteration float64

	for j := 1; j < n; j++ {//sum from j=1 to n-1
		iteration = numbers.MultiplyLog(AlleleFrequencyProbability(j, n, alpha, binomCache), AncestralAscertainmentProbability(n, j, d))
		answer = numbers.AddLog(answer, iteration)
	}

	return answer
}

func DerivedAscertainmentDenominator(alpha float64, n int, d int, binomCache [][]float64) float64 {
	var answer float64 = 0.0
	var iteration float64

	for j := 1; j < n; j++ {//sum from j=1 to n-1
		iteration = numbers.MultiplyLog(AlleleFrequencyProbability(j, n, alpha, binomCache), DerivedAscertainmentProbability(n, j, d))
		answer = numbers.AddLog(answer, iteration)
	}
	return answer
}

func AlleleFrequencyProbabilityDerivedAscertainment(alpha float64, i int, n int, d int, binomCache [][]float64) float64 {
	return numbers.DivideLog(numbers.MultiplyLog(AlleleFrequencyProbability(i, n, alpha, binomCache), DerivedAscertainmentProbability(n, i, d)), DerivedAscertainmentDenominator(alpha, n, d, binomCache))
}

func AlleleFrequencyProbabilityAncestralAscertainment(alpha float64, i int, n int, d int, binomCache [][]float64) float64 {
	return numbers.DivideLog(numbers.MultiplyLog(AlleleFrequencyProbability(i, n, alpha, binomCache), AncestralAscertainmentProbability(n, i, d)), AncestralAscertainmentDenominator(alpha, n, d, binomCache))
}

//AfsLikelihoodDerivedAscertainment is like AFSLikelihood, but makes a correction for divergence-based ascertainment when variant sets were selected for derived alleles between two groups of d individuals.
//More explanation can be found in Katzman et al, this is the inverse of Eq. 11 in the methods.
func AfsLikelihoodDerivedAscertainment(afs AFS, alpha []float64, binomMap [][]float64, d int) float64 {
	var answer float64 = 0.0
	// loop over all segregating sites
	for j := 0; j < len(afs.sites); j++ {
		answer = numbers.MultiplyLog(answer, AlleleFrequencyProbabilityDerivedAscertainment(alpha[j], afs.sites[j].i, afs.sites[j].n, d, binomMap))
	}
	return answer
}

//AfsLikelihoodAncestralAscertainment is like AFSLikelihood, but makes a correction for divergence-based ascertainment when variant sets were selected for ancestral alleles (as in conserved regions like UCEs). d is the number of genomes from each group in the ascertainment process.
func AfsLikelihoodAncestralAscertainment(afs AFS, alpha []float64, binomMap [][]float64, d int) float64 {
	var answer float64 = 0.0
	//loop over all segregating sites
	for j := 0; j < len(afs.sites); j++ {
		answer = numbers.MultiplyLog(answer, AlleleFrequencyProbabilityAncestralAscertainment(alpha[j], afs.sites[j].i, afs.sites[j].n, d, binomMap))
	}
	return answer
}