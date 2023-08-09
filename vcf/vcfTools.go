package vcf

import (
	"strings"
)

// Snp checks SVTYPE and returns true if the value is SNP.
func Snp(v Vcf) bool {
	return strings.Contains(v.Info, "SVTYPE=SNP")
}

// Ins checks SVTYPE and returns true if the value is INS.
func Ins(v Vcf) bool {
	return strings.Contains(v.Info, "SVTYPE=INS")
}

// Del checks SVTYPE and returns true if the value is DEL.
func Del(v Vcf) bool {
	return strings.Contains(v.Info, "SVTYPE=DEL")
}
