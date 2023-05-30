package popgen

import (
	"log"
	"math"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"github.com/vertgenlab/gonomics/vcf"
)

/*
This file implements functions to construct a Hierarchical Bayesian model for inference of the selection parameter
alpha based on an allele frequency spectrum obtained from multiple alignment data.
More details can be obtained from the Ph.D. Thesis of Sol Katzman, available at https://compbio.soe.ucsc.edu/theses/Sol-Katzman-PhD-2010.pdf.
Equations from this thesis that are replicated here are marked.
*/

const IntegralBound float64 = 1e-12

type LikelihoodFunction byte

const (
	Uncorrected LikelihoodFunction = iota
	Ancestral
	Derived
)

// k is len(sites).
type Afs struct {
	Sites []*SegSite
}

// SegSite is the basic struct for segregating sites, used to construct allele frequency spectra.
type SegSite struct {
	I int                //individuals with the allele
	N int                //total number of individuals
	L LikelihoodFunction //specifies with likelihood function to use. 1 for ancestral, 0 for uncorrected, 2 for derived.
}

func segSitesAreEqual(a SegSite, b SegSite) bool {
	if a.I != b.I {
		return false
	}
	if a.N != b.N {
		return false
	}
	if a.L != b.L {
		return false
	}
	return true
}

// InvertSegSite reverses the polarity of a segregating site.
func InvertSegSite(s *SegSite) {
	s.I = s.N - s.I
}

// MultiFaToAfs constructs an allele frequency spectrum struct from a multiFa alignment block.
func MultiFaToAfs(aln []fasta.Fasta) Afs {
	var answer Afs
	var count int
	var current dna.Base
	aln = fasta.SegregatingSites(aln)
	answer.Sites = make([]*SegSite, len(aln[0].Seq))
	for i := 0; i < len(aln[0].Seq); i++ {
		count = 0
		current = aln[0].Seq[i]
		for j := 0; j < len(aln); j++ {
			if aln[j].Seq[i] != current {
				count++
			}
		}
		answer.Sites = append(answer.Sites, &SegSite{count, len(aln), Uncorrected}) //hardcoded to default likelihood for now.
	}
	return answer
}

// VcfToAfs reads in a vcf file, parses the genotype information, and constructs an AFS struct.
// Polarized flag, when true, returns only variants with the ancestor annotated in terms of polarized, derived allele frequencies.
// IncludeRef, when true, uses the reference genome as an additional data point for the output AFS.
func VcfToAfs(filename string, UnPolarized bool, DivergenceAscertainment bool, IncludeRef bool) (*Afs, error) {
	var answer Afs
	answer.Sites = make([]*SegSite, 0)
	alpha, _ := vcf.GoReadToChan(filename)
	var currentSeg *SegSite
	var pass bool
	for i := range alpha {
		currentSeg, pass = VcfSampleToSegSite(i, DivergenceAscertainment, UnPolarized, IncludeRef)
		if pass {
			answer.Sites = append(answer.Sites, currentSeg)
		}
	}
	return &answer, nil
}

// VcfSampleToSegSite returns a SegSite struct from an input Vcf entry. Enables flag for divergenceBasedAscertainment correction conditions and the unPolarized condition.
// Two returns: a pointer to the SegSite struct, and a bool that is true if the SegSite was made without issue, false for soft errors.
// If IncludeRef is true, adds the reference allele as an extra data point to the SegSite.
func VcfSampleToSegSite(i vcf.Vcf, DivergenceAscertainment bool, UnPolarized bool, IncludeRef bool) (*SegSite, bool) {
	var currentSeg = &SegSite{I: 0, N: 0, L: Uncorrected}
	var j int
	//gVCF converts the alt and ref to []DNA.base, so structural variants with <CN0> notation will fail to convert. This check allows us to ignore these cases.
	if !strings.ContainsAny(i.Alt[0], "<>") { //By definition, segregating sites are biallelic, so we only check the first entry in Alt.
		for j = 0; j < len(i.Samples); j++ {
			if len(i.Samples[j].Alleles) == 2 && i.Samples[j].Alleles[0] != -1 && i.Samples[j].Alleles[1] != -1 { //check data for both alleles exist for sample.
				currentSeg.N = currentSeg.N + 2
				if i.Samples[j].Alleles[0] > 0 {
					currentSeg.I++
				}
				if i.Samples[j].Alleles[1] > 0 {
					currentSeg.I++
				}
			}
		}

		if IncludeRef {
			if vcf.IsAltAncestor(i) { //if the reference base is in the derived state.
				currentSeg.I++
			}
			currentSeg.N++
		}

		if currentSeg.N == 0 {
			log.Fatalf("error in VcfToAFS: variant had no sample data")
		}
		if currentSeg.I == 0 || currentSeg.N == currentSeg.I {
			log.Fatalf("error in VcfToAFS: variant is nonsegregating and has an allele frequency of 0 or 1")
		}
		if !UnPolarized && vcf.HasAncestor(i) {
			if vcf.IsRefAncestor(i) && DivergenceAscertainment {
				currentSeg.L = Ancestral
			}
			if vcf.IsAltAncestor(i) {
				InvertSegSite(currentSeg)
				if DivergenceAscertainment {
					currentSeg.L = Derived
				}
			} else if !vcf.IsRefAncestor(i) {
				return nil, false //this special case arises when neither the alt or ref allele is ancestral, can occur with multiallelic positions. In VcfAfs, these positions are not represented in the output.
			}
		}
		if !UnPolarized && !vcf.HasAncestor(i) {
			log.Fatalf("To make a polarized AFS, ancestral alleles must be annotated. Run vcfAncestorAnnotation, filter out variants without ancestral alleles annotated with vcfFilter, or mark unPolarized in options.")
		}
	}

	return currentSeg, true
}

// VcfDerivedAlleleFrequency returns the derived allele frequency based on the sample columns of a VCF variant.
func VcfSampleDerivedAlleleFrequency(v vcf.Vcf) float64 {
	if !vcf.IsPolarizable(v) {
		log.Fatalf("VcfSampleDerivedAlleleFrequency requires polarizable input variants.")
	}
	segSite, _ := VcfSampleToSegSite(v, false, false, false)
	return float64(segSite.I) / float64(segSite.N)
}

// AfsToFrequency converts an  allele frequency spectrum into allele frequencies. Useful for constructing subsequent AFS histograms.
func AfsToFrequency(a Afs) []float64 {
	answer := make([]float64, len(a.Sites))
	for x := 0; x < len(a.Sites); x++ {
		answer[x] = float64(a.Sites[x].I) / float64(a.Sites[x].N)
	}
	return answer
}

// AfsStationarity returns the function value from a stationarity distribution with selection parameter alpha from a particular input allele frequency p.
func AfsStationarity(p float64, alpha float64) float64 {
	return (1 - math.Exp(-alpha*(1-p))) * 2 / ((1 - math.Exp(-alpha)) * p * (1 - p))
}

// AfsStationarityClosure returns a func(float64)float64 for a stationarity distribution with a fixed alpha value for subsequent integration.
func AfsStationarityClosure(alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		return AfsStationarity(p, alpha)
	}
}

// FIntegralComponent is a helper function of AfsSampleDensity and represents the component within the integral.
func FIntegralComponent(n int, k int, alpha float64, binomMap [][]float64) func(float64) float64 {
	var binomCoeff float64 = binomMap[n][k]
	return func(p float64) float64 {
		expression := numbers.BinomialExpressionLog(n-2, k-1, p)
		logPart := math.Log((1 - math.Exp(-alpha*(1.0-p))) * 2 / (1 - math.Exp(-alpha)))
		return logspace.Multiply(binomCoeff, logspace.Multiply(expression, logPart))
	}
}

// AfsSampleDensity (also referred to as the F function) is the product of the stationarity and binomial distributions integrated over p, the allele frequency.
func AfsSampleDensity(n int, k int, alpha float64, binomMap [][]float64, integralError float64) float64 {
	if alpha == 0 {
		log.Fatalf("The stationarity distribution cannot be evaluated with an alpha parameter of exactly zero.")
	}
	var switchPoint float64 = float64(k) / float64(n)
	f := FIntegralComponent(n, k, alpha, binomMap)
	//TODO: Integral accuracy is set at 1e-7, but lowering this may increase runtime without much accuracy cost.
	return logspace.Add(numbers.AdaptiveSimpsonsLog(f, 0.0, switchPoint, integralError, 100), numbers.AdaptiveSimpsonsLog(f, switchPoint, 1.0, integralError, 100))
}

// AlleleFrequencyProbability returns the probability of observing i out of n alleles from a stationarity distribution with selection parameter alpha.
func AlleleFrequencyProbability(i int, n int, alpha float64, binomMap [][]float64, integralError float64) float64 {
	var denominator float64 = math.Inf(-1) //denominator begins at -Inf when in log space
	// j loops over all possible values of i
	for j := 1; j < n; j++ {
		denominator = logspace.Add(denominator, AfsSampleDensity(n, j, alpha, binomMap, integralError))
	}
	return logspace.Divide(AfsSampleDensity(n, i, alpha, binomMap, integralError), denominator)
}

// AfsLikelihood returns P(Data|alpha), or the likelihood of observing a particular allele frequency spectrum given alpha, a vector of selection parameters.
func AfsLikelihood(afs Afs, alpha []float64, binomMap [][]float64, integralError float64) float64 {
	var answer float64 = 0.0
	// loop over all segregating sites
	for j := range afs.Sites {
		answer = logspace.Multiply(answer, AlleleFrequencyProbability(afs.Sites[j].I, afs.Sites[j].N, alpha[j], binomMap, integralError))
	}
	return answer
}

// AfsLikelihoodFixedAlpha calculates the likelihood of observing a particular frequency spectrum for a given alpha, a selection parameter.
// This represents the special case where every segregating site has the same value for selection, which enables faster, more simplified calculation.
func AfsLikelihoodFixedAlpha(afs Afs, alpha float64, binomMap [][]float64, integralError float64) float64 {
	allN := findAllN(afs)
	var answer float64 = 0.0
	likelihoodCache := BuildLikelihoodCache(allN)
	for j := range afs.Sites {
		if likelihoodCache[afs.Sites[j].N][afs.Sites[j].I] == 0.0 { //if this particular segregating site has not already had its likelihood value cached, we want to calculate and cache it.
			likelihoodCache[afs.Sites[j].N][afs.Sites[j].I] = AlleleFrequencyProbability(afs.Sites[j].I, afs.Sites[j].N, alpha, binomMap, integralError)
		}
		answer = logspace.Multiply(answer, likelihoodCache[afs.Sites[j].N][afs.Sites[j].I])
	}
	return answer
}

// AfsLikelihoodFixedAlphaClosure returns a func(float64) float64 representing the likelihood function for a specific derived allele frequency spectrum with a single selection parameter alpha.
func AfsLikelihoodFixedAlphaClosure(afs Afs, binomMap [][]float64, integralError float64) func(float64) float64 {
	return func(alpha float64) float64 {
		return AfsLikelihoodFixedAlpha(afs, alpha, binomMap, integralError)
	}
}

// BuildLikelihoodCache constructs a cache of likelihood values for MLE. likelihoodCache[n][i] stores the likelihood for a segregating site of n alleles with i in the derived state for a particular selection parameter alpha.
func BuildLikelihoodCache(allN []int) [][]float64 {
	answer := make([][]float64, numbers.MaxMany(allN...)+1) //make the first dimension the output matrix large enough to hold the highest observed N.
	for n := range allN {
		answer[allN[n]] = make([]float64, allN[n]) //for each N value in the AFS, set the second dimension to the length of possible i values.
	}
	return answer
}
