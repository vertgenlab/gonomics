package popgen

import (
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"math/rand"
	//DEBUG: "fmt"
)

//SimulateSegSite returns a segregating site with a non-zero allele frequency sampled from a stationarity distribution with selection parameter alpha.
//Divergent flag, which signals a difference between reference genomes, is false by default, but can be inverted to true at later steps.
func SimulateSegSite(alpha float64, n int, boundAlpha float64, boundBeta float64, boundMultiplier float64) *SegSite {
	var fatalCount int = 1000000
	var maxIteration int = 10000000
	var r, x float64
	var count, i int

	bound := numbers.ScaledBetaSampler(boundAlpha, boundBeta, boundMultiplier)
	f := AfsStationarityClosure(alpha)

	for i = 0; i < fatalCount; i++ {
		count = 0
		x, _ = numbers.BoundedRejectionSample(bound, f, 0.0, 1.0, maxIteration)
		for j := 0; j < n; j++ {
			r = rand.Float64()
			if r < x {
				count++
			}
		}
		if count < 1 || count == n {
			continue
		}
		return &SegSite{count, n, false}
	}
	log.Fatalf("Error in simulateSegSite: unable to produce non-zero allele frequency for alpha:%f and %v alleles in 10000 iterations.", alpha, n)
	return &SegSite{0, 0, false}
}

//SimulateGenotype returns a slice of type vcf.GenomeSample, representing a Sample field of a vcf struct, with an allele frequency drawn from a stationarity distribution with selection parameter alpha.
func SimulateGenotype(alpha float64, n int, boundAlpha float64, boundBeta float64, boundMultiplier float64) []vcf.GenomeSample {
	var answer []vcf.GenomeSample = make([]vcf.GenomeSample, 0)
	s := SimulateSegSite(alpha, n, boundAlpha, boundBeta, boundMultiplier)
	alleleArray := SegSiteToAlleleArray(s)
	var d int
	for c := 0; c < n; c += 2 {
		d = c + 1
		//if we have an odd number of alleles, we make one haploid entry
		if d >= n {
			answer = append(answer, vcf.GenomeSample{AlleleOne: alleleArray[c], AlleleTwo: -1, Phased: false, FormatData: make([]string, 0)})
		} else {
			answer = append(answer, vcf.GenomeSample{AlleleOne: alleleArray[c], AlleleTwo: alleleArray[d], Phased: false, FormatData: make([]string, 0)})
		}
	}
	return answer
}

//SegSiteToAlleleArray is a helper function of SimulateGenotype that takes a SegSite, constructs and array of values with i values set to 1 and n-i values set to 0.
//The array is then shuffled and returned.
func SegSiteToAlleleArray(s *SegSite) []int16 {
	var answer []int16 = make([]int16, s.N)
	for j := 0; j < s.I; j++ {
		answer[j] = 1
	}
	rand.Shuffle(len(answer), func(i, j int) { answer[i], answer[j] = answer[j], answer[i] })
	return answer
}

//StationaritySample returns an allele frequency i out of n individuals sampled from a stationarity
//distribution with selection parameter alpha.
func StationaritySampler(alpha float64, samples int, maxSampleDepth int, bins int, xLeft float64, xRight float64) []float64 {
	f := AfsStationarityClosure(alpha)
	return numbers.FastRejectionSampler(xLeft, xRight, f, bins, maxSampleDepth, samples)
}
