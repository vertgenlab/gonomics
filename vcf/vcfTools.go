package vcf

import (
	"strings"
)

func Snp(v Vcf) bool {
	return strings.Contains(v.Info, "SVTYPE=SNP")
}

func Ins(v Vcf) bool {
	return strings.Contains(v.Info, "SVTYPE=INS")
}

func Del(v Vcf) bool {
	return strings.Contains(v.Info, "SVTYPE=DEL")
}
