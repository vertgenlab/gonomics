package popgen

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
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

//AFSStationarityClosure returns a func(float64)float64 for a stationarity distribution with a fixed alpha value for subsequent integration.
func AFSStationarityClosure(alpha float64) func(float64) float64 {
	return func(p float64) float64 {
		return AFSStationarity(p, alpha)
	}
}

/*
//AFSSampleClosure returns a func(float64)float64 for integration based on a stationarity distribution with a fixed alpha selection parameter, sampled with n alleles with k occurances.
func AFSSampleClosure(n int, k int, alpha float64, nkpCache [][][]float64) func(float64) float64 {
	return func(p float64) float64 {
		//fmt.Printf("AFS: %e.\tBinomial:%e\n", AFSStationarity(p, alpha), numbers.BinomialDist(n, k, p))
		return numbers.MultiplyLog(math.Log(AFSStationarity(p, alpha)), binomCache[n][k][p])
	}
}*/

//AFSSAmpleDensity returns the integral of AFSSampleClosure between 0 and 1.
func AFSSampleDensity(n int, k int, alpha float64, nkpCache [][][]float64, alleleFrequencyCache []float64) float64 {
	//first evaluate integral with 1000 bins and 10000 bins and check the error.
	thousandBins := LogIntegrateStationarityCache(alpha, n, k, 100, alleleFrequencyCache, nkpCache)
	tenThousandBins := LogIntegrateStationarityCache(alpha, n, k, 10, alleleFrequencyCache, nkpCache)

	if (thousandBins - tenThousandBins)/tenThousandBins < 1e-8 {
		return tenThousandBins
	}
	return LogIntegrateStationarityCache(alpha, n, k, 1, alleleFrequencyCache, nkpCache)
}

//AlleleFrequencyProbability returns the probability of observing i out of n alleles from a stationarity distribution with selection parameter alpha.
func AlleleFrequencyProbability(i int, n int, alpha float64, nkpCache [][][]float64, alleleFrequencyCache []float64) float64 {
	var denominator float64 = math.Inf(-1)//denominator begins at -Inf when in log space
	//check if n has already been seen
	for j := 1; j < n-1; j++ {
		denominator = numbers.AddLog(denominator, AFSSampleDensity(n, j, alpha, nkpCache, alleleFrequencyCache))
	}
	return numbers.DivideLog(AFSSampleDensity(n, i, alpha, nkpCache, alleleFrequencyCache), denominator)
}

//AfsLikelihood returns P(Data|alpha), or the likelihood of observing a particular allele frequency spectrum given alpha, a vector of selection parameters.
func AFSLikelihood(afs AFS, alpha []float64, nkpCache [][][]float64, alleleFrequencyCache []float64) float64 {
	var answer float64 = 0.0
	for j := 1; j < len(afs.sites); j++ {
		answer = numbers.MultiplyLog(answer, AlleleFrequencyProbability(afs.sites[j].i, afs.sites[j].n, alpha[j], nkpCache, alleleFrequencyCache))
	}
	return answer
}
