package popgen

import (
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"math"
	"strings"
	//DEBUG"fmt"
)

/*
This file implements functions to construct a Hierarchical Bayesian model for inference of the selection parameter
alpha based on an allele frequency spectrum obtained from multiple alignment data.
More details can be obtained from the Ph.D. Thesis of Sol Katzman, available at https://compbio.soe.ucsc.edu/theses/Sol-Katzman-PhD-2010.pdf.
Equations from this thesis that are replicated here are marked.
*/

//k is len(sites)
type AFS struct {
	sites []*SegSite
}

type SegSite struct {
	i int //individuals with the allele
	n int //total number of individuals
}

func GVCFToAFS(filename string) AFS {
	var answer AFS
	answer.sites = make([]*SegSite, 0)
	alpha := vcf.GoReadGVcf(filename)
	var currentSeg *SegSite
	var j int
	for i := range alpha.Vcfs {
		currentSeg = &SegSite{i: 0, n: 0}
		//gVCF converts the alt and ref to []DNA.base, so structural variants with <CN0> notation will fail to convert. This check allows us to ignore these cases.
		if !strings.ContainsAny(i.Alt, "<>") {
			g := vcf.VcfToGvcf(i)
			for j = 0; j < len(g.Genotypes); j++ {
				if g.Genotypes[j].AlleleOne != -1 && g.Genotypes[j].AlleleTwo != -1 { //check data for both alleles exist for sample.
					currentSeg.n = currentSeg.n + 2
					if g.Genotypes[j].AlleleOne > 0 {
						currentSeg.i++
					}
					if g.Genotypes[j].AlleleTwo > 0 {
						currentSeg.i++
					}
				}
			}
			if currentSeg.n != 0 { //catches variants where there is no data from the samples (can happen when filtering columns)
				answer.sites = append(answer.sites, currentSeg)
			}
		}
	}
	alpha.File.Close()
	return answer
}

//converts an  allele frequency spectrum into allele frequencies. Useful for constructing subsequent AFS histograms.
func AFSToFrequency(a AFS) []float64 {
	var answer []float64
	answer = make([]float64, len(a.sites))
	for x := 0; x < len(a.sites); x++ {
		answer[x] = float64(a.sites[x].i) / float64(a.sites[x].n)
	}
	return answer
}

//eq 2.1
func AFSStationarity(p float64, alpha float64) float64 {
	return (1 - math.Exp(-alpha*(1-p))) * 2 / ((1 - math.Exp(-alpha)) * p * (1 - p))
}

func AFSStationarityClosure(alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		return AFSStationarity(p, alpha)
	}
}

func AFSSampleClosure(n int, k int, alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		//fmt.Printf("AFS: %e.\tBinomial:%e\n", AFSStationarity(p, alpha), numbers.BinomialDist(n, k, p))
		return AFSStationarity(p, alpha) * numbers.BinomialDist(n, k, p)
	}
}

//eq. 2.2
func AFSSampleDensity(n int, k int, alpha float64) float64 {
	f := AFSSampleClosure(n, k, alpha)
	//DEBUG prints
	//fmt.Printf("f(0.1)=%e\n", f(0.1))
	//fmt.Printf("AFS: %e.\tBinomial:%e\n", AFSStationarity(0.1, alpha), numbers.BinomialDist(n, k, 0.1))
	//fmt.Printf("N: %v. K: %v. Alpha: %f.\n", n, k, alpha)
	//n choose k * Definiteintegral(p(1-p secrition)stationaritydensity)
	return numbers.DefiniteIntegral(f, 0.000001, 0.9999999999)
}

//eq 2.3
func AlleleFrequencyProbability(i int, n int, alpha float64) float64 {
	var denominator float64
	for j := 1; j < n-1; j++ {
		denominator = denominator + AFSSampleDensity(n, j, alpha)
	}
	return AFSSampleDensity(n, i, alpha) / denominator
}

//eq 2.4
//afs array has a dummy variable in position 0, so loop starts at 1.
func AFSLikelihood(afs AFS, alpha []float64) float64 {
	var answer float64 = 1.0
	for j := 1; j < len(afs.sites); j++ {
		answer = answer * AlleleFrequencyProbability(afs.sites[j].i, afs.sites[j].n, alpha[j])
	}
	return answer
}
