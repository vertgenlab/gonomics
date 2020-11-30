package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
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

func CompareVcf(alpha *Vcf, beta *Vcf) int {
	compareStorage := strings.Compare(alpha.Chr, beta.Chr)
	if compareStorage != 0 {
		return compareStorage
	}
	return CompareCoord(alpha, beta)//TODO: should we also compare genotypes? Would we want to sort more than chr and coord?
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

func CompareGenotypes(alpha []GenomeSample, beta []GenomeSample) int {
	var res int
	stop := numbers.Min(len(alpha), len(beta))
	for i := 0; i < stop; i++ {
		res = CompareGenomeSample(alpha[i], beta[i])
		if res != 0 {
			return res
		}
	}
	if len(alpha) < len(beta) {
		return -1
	} else if len(alpha) > len(beta) {
		return 1
	}
	return 0
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
	if dna.CompareSeqsIgnoreCaseAndGaps(alpha.Ref, beta.Ref) != 0 {
		return false
	}
	if dna.CompareTwoDSeqsIgnoreCaseAndGaps(alpha.Alt, beta.Alt) != 0 {
		return false
	}
	if strings.Compare(alpha.Filter, beta.Filter) != 0 {
		return false
	}
	if strings.Compare(alpha.Info, beta.Info) != 0 {
		return false
	}
	if strings.Compare(alpha.Notes, beta.Notes) != 0 {
		return false
	}
	if CompareGenotypes(alpha.Genotypes, beta.Genotypes) != 0 {
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
