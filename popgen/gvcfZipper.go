package popgen

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"sort"
	"strings"
)

//TajimaZipperVcf returns the Tajima's D value for each entry in a sorted bed slice using a pre-SORTED gVCF file for variant information.
func TajimaZipperVcf(b []*bed.Bed, VcfFile string) []float64 {
	var answer []float64
	var all []*IndividualAllele
	answer = make([]float64, 0)
	alpha, _ := vcf.GoReadToChan(VcfFile) //must be sorted beforehand
	sort.Sort(bed.ByGenomicCoordinates{b})
	var currBedIndex int = 0
	var currVariants []*vcf.Vcf
	currVariants = make([]*vcf.Vcf, 0)
	for i := range alpha {
		for 1 > 0 { //infinite condition for break loop
			if helperVcfOverlapsBed(b[currBedIndex], i) {
				currVariants = append(currVariants, i)
			}
			if i.Pos > int(b[currBedIndex].ChromEnd) {
				all = VariantsToIndividualAlleleSlice(currVariants)
				answer = append(answer, TajimaGVCF(all))
				currVariants = currVariants[:0] //clear the list of current variants
				currBedIndex++
			} else {
				break
			}
		}
	}
	return answer
}

//VariantsToIndividualAlleleSlice converts a slice of vcf records into a slice of IndividualAllele structs.
func VariantsToIndividualAlleleSlice(curr []*vcf.Vcf) []*IndividualAllele {
	var firstTime bool = false
	var j, k int
	var all []*IndividualAllele
	for i := range curr {
		if !strings.ContainsAny(curr[i].Alt[0], "<>") { //gVCF converts the alt and ref to []DNA.base, so structural variants with <CN0> notation will fail to convert. This check allows us to ignore these cases.
			g := vcf.VcfToGvcf(curr[i])
			if firstTime {
				firstTime = false
				all = make([]*IndividualAllele, len(g.Genotypes)*2)
				for k = range all {
					all[k].sites = make([][]dna.Base, 0)
				}
			}
			for j = range g.Genotypes {
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
	return all
}

//helperVcfOverlapsBed checks if a bed and vcf entry are overlapping to avoid imports.
func helperVcfOverlapsBed(b *bed.Bed, v *vcf.Vcf) bool {
	if b.Chrom != v.Chr {
		return false
	}
	return int(b.ChromStart)+1 < v.Pos && int(b.ChromEnd)+1 > v.Pos
}
