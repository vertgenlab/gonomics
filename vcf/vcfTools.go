package vcf

import (
	"strings"
)

func Snp(v Vcf) bool {
	if strings.Contains(v.Info, "SVTYPE=SNP") {
		return true
	}
	return false
}

func Ins(v Vcf) bool {
	if strings.Contains(v.Info, "SVTYPE=INS") {
		return true
	}
	return false
}

func Del(v Vcf) bool {
	var truth bool = false
	if strings.Contains(v.Info, "SVTYPE=DEL") {
		return true
	}
	return truth
}
