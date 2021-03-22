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

//k is len(sites)
type AFS struct {
	sites []*SegSite
}

//SegSite is the basic struct for segregating sites, used to construct allele frequency spectra.
type SegSite struct {
	i int //individuals with the allele
	n int //total number of individuals
}

//InvertSegSite reverses the polarity of a segregating site.
func InvertSegSite(s *SegSite) {
	s.i = s.n - s.i
}

//MultiFaToAFS constructs an allele frequency spectrum struct from a multiFa alignment block.
func MultiFaToAFS(aln []fasta.Fasta) AFS {
	var answer AFS
	var count int
	var current dna.Base
	aln = fasta.SegregatingSites(aln)
	answer.sites = make([]*SegSite, len(aln[0].Seq))
	for i := 0; i < len(aln[0].Seq); i++ {
		count = 0
		current = aln[0].Seq[i]
		for j := 0; j < len(aln); j++ {
			if aln[j].Seq[i] != current {
				count++
			}
		}
		answer.sites = append(answer.sites, &SegSite{count, len(aln)})
	}
	return answer
}

//VcfToAFS reads in a vcf file, parses the genotype information, and constructs an AFS struct.
//Polarized flag, when true, returns only variants with the ancestor annotated in terms of polarized, derived allele frequencies.
func VcfToAFS(filename string, polarized bool) (*AFS, error) {
	var answer AFS
	answer.sites = make([]*SegSite, 0)
	alpha, _ := vcf.GoReadToChan(filename)
	var currentSeg *SegSite
	var j int
	for i := range alpha {
		currentSeg = &SegSite{i: 0, n: 0}
		//gVCF converts the alt and ref to []DNA.base, so structural variants with <CN0> notation will fail to convert. This check allows us to ignore these cases.
		if !strings.ContainsAny(i.Alt[0], "<>") { //By definition, segregating sites are biallelic, so we only check the first entry in Alt.
			for j = 0; j < len(i.Samples); j++ {
				if i.Samples[j].AlleleOne != -1 && i.Samples[j].AlleleTwo != -1 { //check data for both alleles exist for sample.
					currentSeg.n = currentSeg.n + 2
					if i.Samples[j].AlleleOne > 0 {
						currentSeg.i++
					}
					if i.Samples[j].AlleleTwo > 0 {
						currentSeg.i++
					}
				}
			}

			if currentSeg.n == 0 {
				return nil, fmt.Errorf("Error in VcfToAFS: variant had no sample data.")
			}
			if currentSeg.i == 0 || currentSeg.n == currentSeg.i {
				return nil, fmt.Errorf("Error in VcfToAFS: variant is nonsegregating and has an allele frequency of 0 or 1.")
			}
			if polarized && vcf.HasAncestor(i) {
				if vcf.IsAltAncestor(i) {
					InvertSegSite(currentSeg)
				} else if !vcf.IsRefAncestor(i) {
					continue //this special case arises when neither the alt or ref allele is ancestral, can occur with multiallelic positions. For now they are not represented in the output AFS.
				}
			}
			if polarized && !vcf.HasAncestor(i) {
				log.Fatalf("To make a polarized AFS, ancestral alleles must be annotated. Run vcfAncestorAnnotation, filter out variants without ancestral alleles annotated with vcfFilter, or mark unPolarized in options.")
			}
			answer.sites = append(answer.sites, currentSeg)
		}
	}
	return &answer, nil
}

//AFSToFrequency converts an  allele frequency spectrum into allele frequencies. Useful for constructing subsequent AFS histograms.
func AFSToFrequency(a AFS) []float64 {
	var answer []float64
	answer = make([]float64, len(a.sites))
	for x := 0; x < len(a.sites); x++ {
		answer[x] = float64(a.sites[x].i) / float64(a.sites[x].n)
	}
	return answer
}

//AFSStationarity returns the function value from a stationarity distribution with selection parameter alpha from a particular input allele frequency p.
func AFSStationarity(p float64, alpha float64) float64 {
	return (1 - math.Exp(-alpha*(1-p))) * 2 / ((1 - math.Exp(-alpha)) * p * (1 - p))
}

//AFSStationarityClosure returns a func(float64)float64 for a stationarity distribution with a fixed alpha value for subsequent integration.
func AFSStationarityClosure(alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		return AFSStationarity(p, alpha)
	}
}

//FIntegralComponent is a helper function of AFSSampleDensity and represents the component within the integral.
func FIntegralComponent(n int, k int, alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		expression := numbers.BinomialExpressionLog(n-2, k-1, p)
		logPart := math.Log((1 - math.Exp(-alpha*(1.0-p))) * 2 / (1 - math.Exp(-alpha)))
		return numbers.MultiplyLog(expression, logPart)
	}
}

//AFSSampleDensity (also referred to as the F function) is the product of the stationarity and binomial distributions integrated over p, the allele frequency.
func AFSSampleDensity(n int, k int, alpha float64, binomMap [][]float64) float64 {
	var switchPoint float64 = float64(k) / float64(n)
	f := FIntegralComponent(n, k, alpha)
	constantComponent := binomMap[n][k]
	integral := numbers.AddLog(numbers.AdaptiveSimpsonsLog(f, 0.0, switchPoint, 1e-8, 100), numbers.AdaptiveSimpsonsLog(f, switchPoint, 1.0, 1e-8, 100))
	return numbers.MultiplyLog(constantComponent, integral)
}

//AlleleFrequencyProbability returns the probability of observing i out of n alleles from a stationarity distribution with selection parameter alpha.
func AlleleFrequencyProbability(i int, n int, alpha float64, binomMap [][]float64) float64 {
	var denominator float64 = math.Inf(-1) //denominator begins at -Inf when in log space
	// j loops over all possible values of i
	for j := 1; j < n; j++ {
		denominator = numbers.AddLog(denominator, AFSSampleDensity(n, j, alpha, binomMap))
	}
	return numbers.DivideLog(AFSSampleDensity(n, i, alpha, binomMap), denominator)
}

//AfsLikelihood returns P(Data|alpha), or the likelihood of observing a particular allele frequency spectrum given alpha, a vector of selection parameters.
func AFSLikelihood(afs AFS, alpha []float64, binomMap [][]float64) float64 {
	var answer float64 = 0.0
	// loop over all segregating sites
	for j := range afs.sites {
		answer = numbers.MultiplyLog(answer, AlleleFrequencyProbability(afs.sites[j].i, afs.sites[j].n, alpha[j], binomMap))
	}
	return answer
}
