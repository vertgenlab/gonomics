package vcf

import (
	"log"
	"sort"
	"strings"
	"github.com/vertgenlab/gonomics/dna"
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

func EqualGVcf(alpha GVcf, beta GVcf) bool {
	if !isEqual(&alpha.Vcf, &beta.Vcf) {
		return false
	}
	if !EqualGenotypes(alpha.Genotypes, beta.Genotypes) || !EqualSeq(alpha.Seq, beta.Seq) {
		return false
	}
	return true
}

func EqualSeq(alpha [][]dna.Base, beta [][]dna.Base) bool {
	if len(alpha) != len(beta) {
		return false
	}
	for i := 0; i < len(alpha); i++ {
		if dna.CompareSeqsIgnoreCase(alpha[i], beta[i]) != 0 {
			return false
		}
	}
	return true
}

func CompareGenomeSample(alpha GenomeSample, beta GenomeSample) int {
	if alpha.AlleleOne != beta.AlleleOne {
		if alpha.AlleleOne < beta.AlleleOne {
			return -1
		}
		return 1
	}
	if alpha.AlleleTwo != beta.AlleleTwo {
		if alpha.AlleleTwo < beta.AlleleTwo {
			return -1
		}
		return 1
	}
	if alpha.Phased != beta.Phased {
		if !alpha.Phased {
			return -1
		}
		return 1
	}
	return 0
}

func EqualGenotypes(alpha []GenomeSample, beta []GenomeSample) bool {
	if len(alpha) != len(beta) {
		return false
	}
	for i := 0; i < len(alpha); i++ {
		if CompareGenomeSample(alpha[i], beta[i]) != 0 {
			return false
		}
	}
	return true
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
