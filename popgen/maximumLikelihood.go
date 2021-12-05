package popgen

import (
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

//MleSettings delineates the experimental parameters for running maximum likelihood estimation on a set of variants using selectionMLE.
type MleSettings struct {
	Left                    float64 //left bound of MLE search space
	Right                   float64 //right bound of MLE search space
	Error                   float64 //desired accuracy in output maximum likelihood estimate
	UnPolarized             bool    //disable the need for ancestor annotation in the vcf file. Assumes the reference allele is in the ancestral state. Use with caution.
	DivergenceAscertainment bool    //if true, use the divergence-based ascertainment bias corrected likelihood function for selection.
	AncestralDivergenceAscertainment bool //if true, use the divergence-based ascertainment bias corrected likelihood function for selection, but only for bases in the ancestral state.
	DerivedDivergenceAscertainment bool //if true, use the divergence-based ascertainment bias corrected likelihood function for selection, but only for bases in the derived state.
	D                       int     //for DivergenceAscertainment, set the size of the ascertainment subset.
	IntegralError           float64 //set the acceptable error in the internal integral calculations in the likelihood function.
	Verbose                 int     //default 0. When set to 1, debug prints appear in standard output.
}

//SelectionMaximumLikelihoodEstimate performs MLE on an input allele frequency spectrum and writes the result to an output file.
func SelectionMaximumLikelihoodEstimate(data Afs, s MleSettings) float64 {
	allN := findAllN(data)
	binomCache := BuildBinomCache(allN)
	var f func(float64) float64
	if s.AncestralDivergenceAscertainment {
		if s.DerivedDivergenceAscertainment {
			log.Fatalf("Use either the ancestral correction OR the derived correction. If you want to apply both, selected 'DivergenceAscertainment.'")
		} else {
			f = AfsAncestralDivergenceAscertainmentFixedAlphaClosure(data, binomCache, s.D, s.IntegralError)
		}
	} else if s.DerivedDivergenceAscertainment {
		f = AfsDerivedDivergenceAscertainmentFixedAlphaClosure(data, binomCache, s.D, s.IntegralError)
	} else if s.DivergenceAscertainment {
		f = AfsDivergenceAscertainmentFixedAlphaClosure(data, binomCache, s.D, s.IntegralError)
	} else {
		f = AfsLikelihoodFixedAlphaClosure(data, binomCache, s.IntegralError)
	}
	return numbers.GoldenSectionMaxSearch(f, s.Left, s.Right, s.Error)
}
