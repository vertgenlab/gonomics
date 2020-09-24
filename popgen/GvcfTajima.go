package popgen

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	//DEBUG: "fmt"
	"math"
	"strings"
)

/*
This file contains functions for calculating Tajima's D from a gVCF. More information on Tajima's D can be found at https://en.wikipedia.org/wiki/Tajima%27s_D
*/

//IndividualAllele is an internal struct for the Tajima's D calculation from gVCFs. It contains the "column" information of a gVCF block
//In other words, an IndividualAllele represents the state of each segregating site for one allele in a gVCF block.
type IndividualAllele struct {
	sites [][]dna.Base
}

//IndividualDist is a helper function for Tajima's D calculations that calculates the distance (the number of segregating sites that are different by state) between two IndividualAllele structs.
func IndividualDist(a IndividualAllele, b IndividualAllele) int {
	var answer int = 0
	if len(a.sites) != len(b.sites) {
		log.Fatalf("IndividualAllele sites lists must be the same length to calculate distance.")
	}
	for i := 0; i < len(a.sites); i++ {
		if len(a.sites[i]) != len(b.sites[i]) {
			answer++
		} else if dna.Dist(a.sites[i], b.sites[i]) > 0 { //dist requires sequences of equal length, so the first if catches this case.
			answer++
		}
	}
	return answer
}

//GVCFCalculateA2 returns the A2 sum used in the calculation of Tajima's D. More information for this constant can be found at the wikipedia page.
func GVCFCalculateA2(n int) float64 {
	var answer float64

	for i := 1; i < n-1; i++ {
		answer += (1.0 / math.Pow(float64(i), 2))
	}
	return answer
}

//GVCFCalculateA1 returns the A1 sum used in the calculation of Tajima's D. (see wiki)
func GVCFCalculateA1(n int) float64 {
	var answer float64

	for i := 1; i < n-1; i++ {
		answer += (1.0 / float64(i))
	}
	return answer
}

//CalculatePiTajima calculates Pi, the mean pairwise distance between the members of an individual allele matrix.
func CalculatePiTajima(all []*IndividualAllele) float64 {
	var distList []int
	var sum, curr int
	for i := 0; i < len(all); i++ {
		for j := i + 1; j < len(all); j++ { //for all pairs of alleles
			curr = IndividualDist(*all[i], *all[j])
			distList = append(distList, curr)
			sum += curr
		}
	}
	n := len(all)
	for i := 0; i < len(distList); i++ {
		sum += distList[i]
	}
	return (float64(sum) / float64(n))
}

//TajimaGVCF is the main function that takes an individual allele matrix and calculates Tajima's D. Limited to SNPs.
func TajimaGVCF(all []*IndividualAllele) float64 {
	pi := CalculatePiTajima(all)
	n := len(all)
	S := len(all[0].sites)
	a1 := GVCFCalculateA1(n)
	a2 := GVCFCalculateA2(n)
	b1 := (float64(n + 1)) / (3.0 * float64(n-1))
	b2 := 2.0 * (math.Pow(float64(n), 2) + float64(n) + 3.0) / (9.0 * float64(n) * float64(n-1))
	c1 := b1 - (1.0 / a1)
	e1 := c1 / a1
	c2 := b2 - float64(n+2)/(a1*float64(n)) + a2/math.Pow(a1, 2)
	e2 := c2 / (math.Pow(a1, 2) + a2)
	return (pi - (float64(S) / a1)) / math.Sqrt(e1*float64(S)+e2*float64(S)*float64(S-1))
}

//TajimaGVCFBedSet is designed to calculate Tajima's D from a
//gVCF file for variants falling in a particular set of bed regions. Returns one value for the whole set. Limited to SNPs.
func TajimaGVCFBedSet(b []*bed.Bed, gVCFFile string) float64 {
	alpha := vcf.GoReadGVcf(gVCFFile)
	var all []*IndividualAllele
	var j, k int
	var firstTime bool = true
	var intervals []interval.Interval
	intervals = make([]interval.Interval, len(b))
	for i := 0; i < len(b); i++ {
		intervals[i] = b[i]
	}
	tree := interval.BuildTree(intervals)
	for i := range alpha.Vcfs {
		if len(interval.Query(tree, i, "any")) > 0 { //in other words, if the variant overlaps any of the beds
			if !strings.ContainsAny(i.Alt, "<>") { //gVCF converts the alt and ref to []DNA.base, so structural variants with <CN0> notation will fail to convert. This check allows us to ignore these cases.
				g := vcf.VcfToGvcf(i)
				//DEBUG: fmt.Printf("Len: %v.\n", len(g.Genotypes))
				if firstTime {
					firstTime = false
					all = make([]*IndividualAllele, len(g.Genotypes)*2) //makes the individual list, with one entry for each allele

					for k = 0; k < len(all); k++ { //in this loop we initialize all the sites lists
						currSites := make([][]dna.Base, 0)
						all[k] = &IndividualAllele{sites: currSites}
					}
				}
				for j = 0; j < len(g.Genotypes); j++ {
					if g.Genotypes[j].AlleleOne == -1 || g.Genotypes[j].AlleleTwo == -1 { //check that data exists for both alleles
						log.Fatalf("Tajima's D on gVCFs requires complete alignment blocks.")
					} else {
						if g.Genotypes[j].AlleleOne > 0 {
							all[j].sites = append(all[j].sites, g.Seq[1])
						} else {
							all[j].sites = append(all[j].sites, g.Seq[0])
						}
						if g.Genotypes[j].AlleleTwo > 0 {
							all[j+len(g.Genotypes)].sites = append(all[j+len(g.Genotypes)].sites, g.Seq[0])
						} else {
							all[j+len(g.Genotypes)].sites = append(all[j+len(g.Genotypes)].sites, g.Seq[1])
						}
					}
				}
			}
		}
	}
	return TajimaGVCF(all)
}
