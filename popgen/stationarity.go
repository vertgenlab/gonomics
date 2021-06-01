package popgen

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"math"
	"strings"
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
	Uncorrected LikelihoodFunction = 0
	Ancestral   LikelihoodFunction = 1
	Derived     LikelihoodFunction = 2
)

//k is len(sites)
type Afs struct {
	Sites []*SegSite
}

//SegSite is the basic struct for segregating sites, used to construct allele frequency spectra.
type SegSite struct {
	I int                //individuals with the allele
	N int                //total number of individuals
	L LikelihoodFunction //specifies with likelihood function to use. 1 for ancestral, 0 for uncorrected, 2 for derived.
}

//InvertSegSite reverses the polarity of a segregating site.
func InvertSegSite(s *SegSite) {
	s.I = s.N - s.I
}

//MultiFaToAfs constructs an allele frequency spectrum struct from a multiFa alignment block.
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

//VcfToAfs reads in a vcf file, parses the genotype information, and constructs an AFS struct.
//Polarized flag, when true, returns only variants with the ancestor annotated in terms of polarized, derived allele frequencies.
func VcfToAfs(filename string, UnPolarized bool, DivergenceAscertainment bool) (*Afs, error) {
	var answer Afs
	answer.Sites = make([]*SegSite, 0)
	alpha, _ := vcf.GoReadToChan(filename)
	var currentSeg *SegSite
	var j int
	for i := range alpha {
		currentSeg = &SegSite{I: 0, N: 0, L: Uncorrected}
		//gVCF converts the alt and ref to []DNA.base, so structural variants with <CN0> notation will fail to convert. This check allows us to ignore these cases.
		if !strings.ContainsAny(i.Alt[0], "<>") { //By definition, segregating sites are biallelic, so we only check the first entry in Alt.
			for j = 0; j < len(i.Samples); j++ {
				if i.Samples[j].AlleleOne != -1 && i.Samples[j].AlleleTwo != -1 { //check data for both alleles exist for sample.
					currentSeg.N = currentSeg.N + 2
					if i.Samples[j].AlleleOne > 0 {
						currentSeg.I++
					}
					if i.Samples[j].AlleleTwo > 0 {
						currentSeg.I++
					}
				}
			}

			if currentSeg.N == 0 {
				return nil, fmt.Errorf("error in VcfToAFS: variant had no sample data")
			}
			if currentSeg.I == 0 || currentSeg.N == currentSeg.I {
				return nil, fmt.Errorf("error in VcfToAFS: variant is nonsegregating and has an allele frequency of 0 or 1")
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
					continue //this special case arises when neither the alt or ref allele is ancestral, can occur with multiallelic positions. For now they are not represented in the output AFS.
				}
			}
			if !UnPolarized && !vcf.HasAncestor(i) {
				log.Fatalf("To make a polarized AFS, ancestral alleles must be annotated. Run vcfAncestorAnnotation, filter out variants without ancestral alleles annotated with vcfFilter, or mark unPolarized in options.")
			}
			answer.Sites = append(answer.Sites, currentSeg)
		}
	}
	return &answer, nil
}

//AfsToFrequency converts an  allele frequency spectrum into allele frequencies. Useful for constructing subsequent AFS histograms.
func AfsToFrequency(a Afs) []float64 {
	var answer []float64
	answer = make([]float64, len(a.Sites))
	for x := 0; x < len(a.Sites); x++ {
		answer[x] = float64(a.Sites[x].I) / float64(a.Sites[x].N)
	}
	return answer
}

//AfsStationarity returns the function value from a stationarity distribution with selection parameter alpha from a particular input allele frequency p.
func AfsStationarity(p float64, alpha float64) float64 {
	return (1 - math.Exp(-alpha*(1-p))) * 2 / ((1 - math.Exp(-alpha)) * p * (1 - p))
}

//AfsStationarityClosure returns a func(float64)float64 for a stationarity distribution with a fixed alpha value for subsequent integration.
func AfsStationarityClosure(alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		return AfsStationarity(p, alpha)
	}
}

//FIntegralComponent is a helper function of AfsSampleDensity and represents the component within the integral.
func FIntegralComponent(n int, k int, alpha float64, binomMap [][]float64) func(float64) float64 {
	var binomCoeff float64 = binomMap[n][k]
	return func(p float64) float64 {
		expression := numbers.BinomialExpressionLog(n-2, k-1, p)
		logPart := math.Log((1 - math.Exp(-alpha*(1.0-p))) * 2 / (1 - math.Exp(-alpha)))
		return numbers.MultiplyLog(binomCoeff, numbers.MultiplyLog(expression, logPart))
	}
}

//AfsSampleDensity (also referred to as the F function) is the product of the stationarity and binomial distributions integrated over p, the allele frequency.
func AfsSampleDensity(n int, k int, alpha float64, binomMap [][]float64, integralError float64) float64 {
	var switchPoint float64 = float64(k) / float64(n)
	f := FIntegralComponent(n, k, alpha, binomMap)
	//TODO: Integral accuracy is set at 1e-7, but lowering this may increase runtime without much accuracy cost.
	return numbers.AddLog(numbers.AdaptiveSimpsonsLog(f, 0.0, switchPoint, integralError, 100), numbers.AdaptiveSimpsonsLog(f, switchPoint, 1.0, integralError, 100))
}

//AlleleFrequencyProbability returns the probability of observing i out of n alleles from a stationarity distribution with selection parameter alpha.
func AlleleFrequencyProbability(i int, n int, alpha float64, binomMap [][]float64, integralError float64) float64 {
	var denominator float64 = math.Inf(-1) //denominator begins at -Inf when in log space
	// j loops over all possible values of i
	for j := 1; j < n; j++ {
		denominator = numbers.AddLog(denominator, AfsSampleDensity(n, j, alpha, binomMap, integralError))
	}
	return numbers.DivideLog(AfsSampleDensity(n, i, alpha, binomMap, integralError), denominator)
}

//AfsLikelihood returns P(Data|alpha), or the likelihood of observing a particular allele frequency spectrum given alpha, a vector of selection parameters.
func AfsLikelihood(afs Afs, alpha []float64, binomMap [][]float64, integralError float64) float64 {
	var answer float64 = 0.0
	// loop over all segregating sites
	for j := range afs.Sites {
		answer = numbers.MultiplyLog(answer, AlleleFrequencyProbability(afs.Sites[j].I, afs.Sites[j].N, alpha[j], binomMap, integralError))
	}
	return answer
}
