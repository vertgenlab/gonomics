package gtf

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	//"github.com/vertgenlab/gonomics/dna"
	//"github.com/vertgenlab/gonomics/fasta".
	"github.com/vertgenlab/gonomics/vcf"
)

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

func FilterVariantExon(v *vcf.Vcf, g map[string]*Gene, c map[string]*chromInfo.ChromInfo) bool {
	a := ExonBoolArray(g, c)
	return VariantArrayOverlap(v, a)
}

func FilterVariantCds(v *vcf.Vcf, g map[string]*Gene, c map[string]*chromInfo.ChromInfo) bool {
	a := CdsBoolArray(g, c)
	return VariantArrayOverlap(v, a)
}

func FilterVariantThreeUtr(v *vcf.Vcf, g map[string]*Gene, c map[string]*chromInfo.ChromInfo) bool {
	a := ThreeUtrBoolArray(g, c)
	return VariantArrayOverlap(v, a)
}

func FilterVariantFiveUtr(v *vcf.Vcf, g map[string]*Gene, c map[string]*chromInfo.ChromInfo) bool {
	a := FiveUtrBoolArray(g, c)
	return VariantArrayOverlap(v, a)
}
