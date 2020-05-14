package vcf

import (
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

func sameRecord(a *Vcf, b *Vcf) bool {
	if isEqual(a, b) {
		return true
	}
	if strings.Compare(a.Chr, b.Chr) == 0 && a.Pos == b.Pos {
		if strings.Compare(a.Ref, b.Ref) == 0 && strings.Compare(a.Alt, b.Alt) == 0 {
			return true
		}
	}
	return false
}
