package gtf

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
)

// FilterVariantGtf takes a vcf record, a gene list from a gtf, ChromInfo to know the length of chromosomes, and whether the function
// should search for exon, coding, 5' UTR, or 3' UTR overlaps.  If more than one type of overlap is selected by setting the multiple of:
// exon, cds, five, three to true, the function returns the logical or of whether the vcf record overlaps that function.
func FilterVariantGtf(v *vcf.Vcf, g map[string]*Gene, c map[string]*chromInfo.ChromInfo, exon bool, code bool, five bool, three bool) bool {
	if exon {
		if !FilterVariantExon(v, g, c) {
			return false
		}
	}
	if code {
		if !FilterVariantCds(v, g, c) {
			return false
		}
	}
	if three {
		if !FilterVariantThreeUtr(v, g, c) {
			return false
		}
	}
	if five {
		if !FilterVariantFiveUtr(v, g, c) {
			return false
		}
	}
	return true
}

// FilterVariantExon take a vcf record, a gene list from a gtf, and ChromInfo to know the length of chromosomes.
// The function returns true if the vcf record overlaps an exon in the gtf.
func FilterVariantExon(v *vcf.Vcf, g map[string]*Gene, c map[string]*chromInfo.ChromInfo) bool {
	a := ExonBoolArray(g, c)
	return VariantArrayOverlap(v, a)
}

// FilterVariantCds take a vcf record, a gene list from a gtf, and ChromInfo to know the length of chromosomes.
// The function returns true if the vcf record overlaps a cds (protein-coding sequence) in the gtf.
func FilterVariantCds(v *vcf.Vcf, g map[string]*Gene, c map[string]*chromInfo.ChromInfo) bool {
	a := CdsBoolArray(g, c)
	return VariantArrayOverlap(v, a)
}

// FilterVariantThreeUtr take a vcf record, a gene list from a gtf, and ChromInfo to know the length of chromosomes.
// The function returns true if the vcf record overlaps a 3' UTR in the gtf.
func FilterVariantThreeUtr(v *vcf.Vcf, g map[string]*Gene, c map[string]*chromInfo.ChromInfo) bool {
	a := ThreeUtrBoolArray(g, c)
	return VariantArrayOverlap(v, a)
}

// FilterVariantFiveUtr take a vcf record, a gene list from a gtf, and ChromInfo to know the length of chromosomes.
// The function returns true if the vcf record overlaps a 5' UTR in the gtf.
func FilterVariantFiveUtr(v *vcf.Vcf, g map[string]*Gene, c map[string]*chromInfo.ChromInfo) bool {
	a := FiveUtrBoolArray(g, c)
	return VariantArrayOverlap(v, a)
}

func FindPromoter(genes []string, upstream int, downstream int, gtf map[string]*Gene, size map[string]chromInfo.ChromInfo) []bed.Bed {
	var answer []bed.Bed
	var currGene, transcript int
	var exists bool
	var info *Gene
	var name string
	var trans *Transcript
	var newBed = bed.Bed{Chrom: "", ChromStart: 0, ChromEnd: 0, FieldsInitialized: 4}

	for currGene = range genes {
		name = genes[currGene]

		info, exists = gtf[name]
		if exists {
			for transcript = range info.Transcripts {
				trans = info.Transcripts[transcript]
				if trans.Strand {
					newBed = bed.Bed{Chrom: trans.Chr, ChromStart: numbers.Max(trans.Start-upstream, 0), ChromEnd: numbers.Min(trans.Start+downstream, size[trans.Chr].Size), Name: info.GeneName, FieldsInitialized: 4}
				} else if !trans.Strand {
					newBed = bed.Bed{Chrom: trans.Chr, ChromStart: trans.Start + downstream, ChromEnd: trans.Start + upstream}
				}
				answer = append(answer, newBed)
			}
		}
	}

	return answer
}
