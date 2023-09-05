package gtf

import (
	"github.com/vertgenlab/gonomics/chromInfo"
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
