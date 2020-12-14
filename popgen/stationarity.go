package popgen

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"math"
	//DEBUG: "log"
	"strings"
	//DEBUG"fmt"
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

type SegSite struct {
	i int //individuals with the allele
	n int //total number of individuals
}

//MultiFaToAFS constructs an allele frequency spectrum struct from a multiFa alignment block.
//TODO: Ask Craig about derived state here.
func MultiFaToAFS(aln []*fasta.Fasta) AFS {
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

//GvcfToAFS reads in a Gvcf file, parses the genotype information, and constructs an AFS struct.
//TODO: This function will change when we update the gVCF stuff.
func VcfToAFS(filename string) AFS {
	var answer AFS
	answer.sites = make([]*SegSite, 0)
	alpha, _ := vcf.GoReadToChan(filename)
	var currentSeg *SegSite
	var j int
	for i := range alpha {
		currentSeg = &SegSite{i: 0, n: 0}
		//gVCF converts the alt and ref to []DNA.base, so structural variants with <CN0> notation will fail to convert. This check allows us to ignore these cases.
		if !strings.ContainsAny(i.Alt[0], "<>") {//By definition, segregting sites are biallelic, so we only check the first entry in Alt.
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
			if currentSeg.n != 0 { //catches variants where there is no data from the samples (can happen when filtering columns)
				answer.sites = append(answer.sites, currentSeg)
			}
		}
	}
	return answer
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

func DetectionProbability(p float64, n int) float64 {
	var pNotDetected float64 = numbers.AddLog(numbers.BinomialExpressionLog(n, 0, p), numbers.BinomialExpressionLog(n, n, p))
	//DEBUG: log.Printf("pNotDetected: %f.", pNotDetected)
	return numbers.SubtractLog(0, pNotDetected)
}

func AfsStationarityCorrected(p float64, alpha float64, n int) float64 {
	//DEBUG: log.Printf("p: %f. n: %d. alpha: %f. Detection: %f. Stationarity: %f.", p, n, alpha, math.Exp(DetectionProbability(p, n)), math.Exp(AFSStationarity(p, alpha)))
	return numbers.MultiplyLog(DetectionProbability(p, n), AFSStationarity(p, alpha))
}

//AFSStationarityClosure returns a func(float64)float64 for a stationarity distribution with a fixed alpha value for subsequent integration.
func AFSStationarityClosure(alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		return AFSStationarity(p, alpha)
	}
}

func AfsSampleClosure(n int, k int, alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		return numbers.MultiplyLog(math.Log(AFSStationarity(p, alpha)), numbers.BinomialExpressionLog(n, k, p))
	}
}

func FIntegralComponent(n int, k int, alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		expression := numbers.BinomialExpressionLog(n-2, k-1, p)
		logPart := math.Log((1-math.Exp(-alpha*(1.0-p))) * 2 / (1-math.Exp(-alpha)))
		//log.Printf("Expression: %v. LogPart: %v.", expression, logPart)
		return numbers.MultiplyLog(expression, logPart)
	}
}

//AFSSampleClosure returns a func(float64)float64 for integration based on a stationarity distribution with a fixed alpha selection parameter, sampled with n alleles with k occurances.
func AFSSampleClosureOld(n int, k int, alpha float64, binomMap [][]float64) func(float64) float64 {
	return func(p float64) float64 {
		//DEBUG: fmt.Println(binomMap)
		//fmt.Printf("AFS: %e.\tBinomial:%e\n", AFSStationarity(p, alpha), numbers.BinomialDist(n, k, p))
		return numbers.MultiplyLog(math.Log(AFSStationarity(p, alpha)), numbers.BinomialDistLogSlice(n, k, p, binomMap))
	}
}

func AFSSampleDensity(n int, k int, alpha float64, binomMap [][]float64) float64 {
	//DEBUG: log.Printf("n: %d. k: %d. alpha: %v.", n, k, alpha)
	f := FIntegralComponent(n, k, alpha)
	//log.Printf("f(0): %f. f(0.25): %f. f(0.5): %f. f(1): %f.", f(0.0), f(0.25), f(0.5), f(1))
	//log.Fatal()
	//constantComponent := numbers.MultiplyLog(binomMap[n][k], math.Log(2 / (1-math.Exp(-alpha))))
	constantComponent := binomMap[n][k]
	return numbers.MultiplyLog(constantComponent, numbers.AdaptiveSimpsonsLog(f, 0.0, 1.0, 1e-8, 100))
}

//AFSSAmpleDensity returns the integral of AFSSampleClosure between 0 and 1.
func AFSSampleDensityOld(n int, k int, alpha float64, binomMap [][]float64) float64 {
	f := AFSSampleClosureOld(n, k, alpha, binomMap)
	//DEBUG prints
	//fmt.Printf("f(0.1)=%e\n", f(0.1))
	//fmt.Printf("AFS: %e.\tBinomial:%e\n", AFSStationarity(0.1, alpha), numbers.BinomialDist(n, k, 0.1))
	//fmt.Printf("N: %v. K: %v. Alpha: %f.\n", n, k, alpha)
	//n choose k * Definiteintegral(p(1-p secrition)stationaritydensity)
	return numbers.LogIntegrateIterative(f, 0.000001, 0.9999999, 20, 10e-8)
}

//AlleleFrequencyProbability returns the probability of observing i out of n alleles from a stationarity distribution with selection parameter alpha.
func AlleleFrequencyProbability(i int, n int, alpha float64, binomMap [][]float64) float64 {
	var denominator float64 = math.Inf(-1)//denominator begins at -Inf when in log space
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
	for j := 0; j < len(afs.sites); j++ {
		answer = numbers.MultiplyLog(answer, AlleleFrequencyProbability(afs.sites[j].i, afs.sites[j].n, alpha[j], binomMap))
	}
	return answer
}
