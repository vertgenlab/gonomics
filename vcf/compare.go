package vcf

import (
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
