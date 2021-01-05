package vcf

import (
	//"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"sort"
	"strings"
)

//CompareCoord compares two VCF structs by Pos for sorting or equality testing.
func CompareCoord(alpha *Vcf, beta *Vcf) int {
	if alpha.Pos < beta.Pos {
		return -1
	}
	if alpha.Pos > beta.Pos {
		return 1
	}
	return 0
}

//CompareVcf compares two Vcf structs for sorting or equality testing.
func CompareVcf(alpha *Vcf, beta *Vcf) int {
	compareStorage := strings.Compare(alpha.Chr, beta.Chr)
	if compareStorage != 0 {
		return compareStorage
	}
	return CompareCoord(alpha, beta) //TODO: should we also compare genotypes? Would we want to sort more than chr and coord?
}

//CompareGenomeSample compares two GenomeSample structs for sorting or equality testing.
func CompareGenomeSample(alpha GenomeSample, beta GenomeSample) int {
	var res int
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
	res = CompareFormatData(alpha.FormatData, beta.FormatData)
	if res != 0 {
		return res
	}
	return 0
}

//CompareFormatData compares the FormatData field of a VCF GenomeSample
func CompareFormatData(alpha []string, beta []string) int {
	return CompareAlt(alpha, beta) //both functions compare slice of strings for equality, but I made this a separate function for readability
}

//CompareGenoeSamples compares a slice of GenomeSample structs which underlies the VCF struct
func CompareGenomeSamples(alpha []GenomeSample, beta []GenomeSample) int {
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

//CompareAlt compares the two slice of string Alt fields from a VCF lexigraphically.
func CompareAlt(alpha []string, beta []string) int {
	var res int
	stop := numbers.Min(len(alpha), len(beta))
	for i := 0; i < stop; i++ {
		res = strings.Compare(alpha[i], beta[i])
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

//Sort sorts a slice of Vcf structs in place.
func Sort(vcfFile []*Vcf) {
	sort.Slice(vcfFile, func(i, j int) bool { return CompareVcf(vcfFile[i], vcfFile[j]) == -1 })
}

//isEqual returns true if two input Vcf structs contain identical information, false otherwise.
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
	if CompareAlt(alpha.Alt, beta.Alt) != 0 {
		return false
	}
	if strings.Compare(alpha.Filter, beta.Filter) != 0 {
		return false
	}
	if strings.Compare(alpha.Info, beta.Info) != 0 {
		return false
	}
	if CompareGenomeSamples(alpha.Samples, beta.Samples) != 0 {
		return false
	}
	return true

}

//AllEqual returns true if each Vcf in a slice of Vcf structs contain identical information.
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

//LogEqual prints which two corresponding entries in two slices of Vcf structs are unequal and fatals out.
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
