package popgen

import (
	"sort"
	"strings"
	"log"
	"github.com/vertgenlab/gonomics/vcf"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
)

//TajimaZipperVcf returns the Tajima's D value for each entry in a sorted bed slice using a pre-SORTED gVCF file for variant information.
func TajimaZipperVcf(b []*bed.Bed, gVCFFile string) []float64 {
	var answer []float64
	var all []*IndividualAllele
	answer = make([]float64, 0)
	alpha := vcf.GoReadGVcf(gVCFFile)//must be sorted beforehand
	sort.Sort(bed.ByGenomicCoordinates{b})
	var currBedIndex int = 0
	var currVariants []*vcf.Vcf
	currVariants = make([]*vcf.Vcf, 0)
	for i := range alpha.Vcfs {
		for 1 > 0 {//infinite condition for break loop
			if helperVcfOverlapsBed(b[currBedIndex], i) {
				currVariants = append(currVariants, i)
			}
			if i.Pos > b[currBedIndex].ChromEnd {
				all = VariantsToIndividualAlleleSlice(currVariants)
				answer = append(answer, TajimaGVCF(all))
				currVariants = currVariants[:0]//clear the list of current variants
				currBedIndex++
			} else {
				break
			}
		}
	}
	return answer
}

func VariantsToIndividualAlleleSlice(curr []*vcf.Vcf) []*IndividualAllele {
	var firstTime bool = false
	var j int
	var all []*IndividualAllele
	for i := 0; i < len(curr); i++ {
		if !strings.ContainsAny(curr[i].Alt, "<>") { //gVCF converts the alt and ref to []DNA.base, so structural variants with <CN0> notation will fail to convert. This check allows us to ignore these cases.
			g := vcf.VcfToGvcf(curr[i])
			if firstTime {
				firstTime = false
				all = make([]*IndividualAllele, len(g.Genotypes)*2)
				for k := 0; k < len(all); k++ {
					all[k].sites = make([][]dna.Base, 0)
				}
			}
			for j = 0; j < len(g.Genotypes); j++ {
				if g.Genotypes[j].AlleleOne == -1 || g.Genotypes[j].AlleleTwo == -1 {//check that data exists for both alleles
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

func helperVcfOverlapsBed(b *bed.Bed, v *vcf.Vcf) bool {
	if b.Chrom != v.Chr {
		return false
	}
	return b.ChromStart + 1 < v.Pos && b.ChromEnd + 1 > v.Pos
}
