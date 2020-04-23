package vcf

import (
	"log"
	"sort"
	"strings"
)

func CompareCoord(alpha *Vcf, beta *Vcf) int {
	if alpha.Pos < beta.Pos {
		return -1
	}
	if alpha.Pos > beta.Pos {
		return 1
	}
	return 0
}

func CompareName(alpha string, beta string) int {
	return strings.Compare(alpha, beta)
}

func CompareVcf(alpha *Vcf, beta *Vcf) int {
	compareStorage := CompareName(alpha.Chr, beta.Chr)
	if compareStorage != 0 {
		return compareStorage
	} else {
		return CompareCoord(alpha, beta)
	}
}

func Sort(vcfFile []*Vcf) {
	sort.Slice(vcfFile, func(i, j int) bool { return CompareVcf(vcfFile[i], vcfFile[j]) == -1 })
}

func isEqual(alpha *Vcf, beta *Vcf) bool {
	if strings.Compare(alpha.Chr, beta.Chr) != 0 {
		return false
	}
	if alpha.Pos != beta.Pos {
		return false
	}
	if strings.Compare(alpha.Id, beta.Id) != 0 {
		return false
	}
	if strings.Compare(alpha.Ref, beta.Ref) != 0 {
		return false
	}
	if strings.Compare(alpha.Alt, beta.Alt) != 0 {
		return false
	}
	if strings.Compare(alpha.Filter, beta.Filter) != 0 {
		return false
	}
	if strings.Compare(alpha.Info, beta.Info) != 0 {
		return false
	}
	return true

}

func AllEqual(alpha []*Vcf, beta []*Vcf) bool {
	if len(alpha) != len(beta) {
		return false
	}
	for i := 0; i < len(alpha); i++ {
		if !isEqual(alpha[i], beta[i]) {
			return false
		}
	}
	return true
}

func LogEqual(alpha []*Vcf, beta []*Vcf) {
	if len(alpha) != len(beta) {
		log.Fatalf("len=%v and len=%v are not equal\n", len(alpha), len(beta))
	}
	for i := 0; i < len(alpha); i++ {
		if !isEqual(alpha[i], beta[i]) {
			log.Fatalf("%v and %v are not equal\n", alpha[i], beta[i])
		}
	}
}

func LowQual(hap Haplotype) bool {
	if hap.One < 0 || hap.Two < 0 {
		return true
	} else {
		return false
	}
}

func Unqualified(key []int16, gt []Haplotype) bool {
	for i := 0; i < len(key); i++ {
		if LowQual(gt[key[i]]) {
			return true
		}
	}
	return false
}

func IsHeterozygous(hap Haplotype) bool {
	if hap.One != hap.Two {
		return true
	}
	if hap.One == hap.Two {
		return false
	}
	return false
}

func IsHomozygous(hap Haplotype) bool {
	if hap.One == hap.Two {
		return true
	}
	if hap.One != hap.Two {
		return false
	}
	
	return false
}

func Heterozygous(key []int16, gt []Haplotype) bool {
	for _, Aa := range key {
		if !IsHeterozygous(gt[Aa]) {
			return false
		}
	}
	return true
}
//Both Homozygous, but different allele
func UniqueHomozygous(key []int16, AA []Haplotype) bool {
	hash := make(map[int16]bool)
	var ok bool
	for aa := 0; aa < len(key); aa++ {
		if IsHomozygous(AA[key[aa]]) {
			_, ok = hash[AA[key[aa]].One]
			if !ok {
				hash[AA[key[aa]].One] = true
			} else {
				return false
			}	
		} else {
			return false
		}
	}
	return true
}

func homozygousPairDiff(AA Haplotype, aa Haplotype) bool {
	if AA.One != aa.Two && IsHomozygous(AA) && IsHomozygous(aa) {
		return true
	} else {
		return false
	}
}