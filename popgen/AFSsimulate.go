package popgen

import (
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"math/rand"
	//DEBUG: "fmt"
)

//SimulateAFS returns an allele frequency spectrum AFS struct for n haploid individuals with k segregating sites
//with a selection parameter alpha.
func SimulateAFS(alpha float64, n int, k int) AFS {
	var answer AFS
	answer.sites = make([]*SegSite, 0)
	var alleleFrequencies []float64 = StationaritySampler(alpha, k, 100000, 1000, 0.001, 0.999, true)
	var count int
	var r float64
	for i := 0; i < k; i++ {
		//now we simulate the discrete number of alleles for n individuals based on the allele frequency
		count = 0//this variable represents the number of times a derived allele is observed
		for j := 0; j < n; j++ {
			r = rand.Float64()
			if r < alleleFrequencies[i] {
				count++ 
			}
		}
		answer.sites = append(answer.sites, &SegSite{count, n})
	}
	return answer
}

//SimulateGenotype creates a matrix of GenomeSamples with k rows and n/2 columns, corresponding to k segregating sites and n alleles.
//based on a selection parameter alpha. Alleles are shuffled among the individuals.
func SimulateGenotype(alpha float64, n int, k int) [][]vcf.GenomeSample {
	a := SimulateAFS(alpha, n, k)
	var answer [][]vcf.GenomeSample = make([][]vcf.GenomeSample, 0, k)
	
	//this memory allocation is no longer needed
	/*for i := 0; i < len(answer); i++ {
		answer[i] = make([]vcf.GenomeSample, 0, len(a.sites)/2 + 1)
	}*/

	var alleleArray []int16
	var d int
	var currRow []vcf.GenomeSample = make([]vcf.GenomeSample, 0, k)
	//var currRow []vcf.GenomeSample
	//for each segregating site, we fill in the  matrix (dimensions r by c)
	for r := 0; r < k; r++ {
		alleleArray = SegSiteToAlleleArray(a.sites[r])
		currRow = make([]vcf.GenomeSample, 0, k)

		for c := 0; c < n; c += 2 {
			d = c + 1
			//if we have an odd number of alleles, we make one haploid entry
			if d >= n {
					currRow = append(currRow, vcf.GenomeSample{AlleleOne:alleleArray[c], AlleleTwo: -1, Phased: false})
				} else {
					currRow = append(currRow, vcf.GenomeSample{AlleleOne: alleleArray[c], AlleleTwo: alleleArray[d], Phased: false})
				}
			//DEBUG: fmt.Printf("%v\n", currRow)
		}
		answer = append(answer, currRow)
	}
	/*DEBUG:
	for i := 0; i < len(answer); i++ {
		fmt.Printf("%v\n", answer[i])
	}*/

	return answer
}

//SegSiteToAlleleArray is a helper function of SimulateGenotype that takes a SegSite, constructs and array of values with i values set to 1 and n-i values set to 0.
//The array is then shuffled and returned.
func SegSiteToAlleleArray(s *SegSite) []int16 {
	var answer []int16 = make([]int16, s.n)
	for j := 0; j < s.i; j++ {
		answer[j] = 1
	}
	rand.Shuffle(len(answer), func(i, j int) { answer[i], answer[j] = answer[j], answer[i]})
	return answer
} 

/*
func AfsToGenotype(a AFS, phase bool) []vcf.GenomeSample {
	var answer []vcf.GenomeSample = make([]vcf.GenomeSample, len(a.sites)/2 + 1)
	for i := 0; i < len(a.sites); i++ {

	}
	return answer
}*/

//StationaritySample returns an allele frequency i out of n individuals sampled from a stationarity 
//distribution with selection parameter alpha.
func StationaritySampler(alpha float64, samples int, maxSampleDepth int, bins int, xLeft float64, xRight float64, randSeed bool) []float64 {
	f := AFSStationarityClosure(alpha)
	return numbers.FastRejectionSampler(xLeft, xRight, f, bins, maxSampleDepth, samples, randSeed)
}
