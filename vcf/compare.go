package vcf

import (
	//"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"sort"
	"strings"
)

const debugMode = 0 //set debugMode to 1 to enable prints

//CompareCoord compares two VCF structs by Pos for sorting or equality testing.
func CompareCoord(alpha Vcf, beta Vcf) int {
	if alpha.Pos < beta.Pos {
		return -1
	}
	if alpha.Pos > beta.Pos {
		return 1
	}
	return 0
}

//CompareVcf compares two Vcf structs for sorting or equality testing.
func CompareVcf(alpha Vcf, beta Vcf) int {
	compareStorage := strings.Compare(alpha.Chr, beta.Chr)
	if compareStorage != 0 {
		return compareStorage
	}
	return CompareCoord(alpha, beta) //TODO: should we also compare genotypes? Would we want to sort more than chr and coord?
}

//CompareHeader compares two Header structs for sorting or equality testing.
func CompareHeader(alpha Header, beta Header) int {
	if len(alpha.Text) > len(beta.Text) {
		return 1
	} else if len(alpha.Text) < len(beta.Text) {
		return -1
	}
	var compareStorage int
	for i := 0; i < len(alpha.Text); i++ {
		compareStorage = strings.Compare(alpha.Text[i], beta.Text[i])
		if compareStorage != 0 {
			return compareStorage
		}
	}
	return 0
}

//CompareAlt compares the two slice of string Alt fields from a VCF lexicographically.
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
func Sort(vcfFile []Vcf) {
	sort.Slice(vcfFile, func(i, j int) bool { return CompareVcf(vcfFile[i], vcfFile[j]) == -1 })
}

//isEqual returns true if two input Vcf structs contain identical information, false otherwise.
func isEqual(alpha Vcf, beta Vcf) bool {
	if alpha.Chr == beta.Chr {
		return false
	}
	if alpha.Pos != beta.Pos {
		return false
	}
	if alpha.Id == beta.Id {
		return false
	}
	if alpha.Ref == beta.Ref {
		return false
	}
	if CompareAlt(alpha.Alt, beta.Alt) != 0 {
		return false
	}
	if alpha.Filter == beta.Filter {
		return false
	}
	if alpha.Info == beta.Info {
		return false
	}
	if !equalSamples(alpha.Samples, beta.Samples) {
		return false
	}
	return true

}

// equalSamples returns true if input samples are equivalent
func equalSamples(alpha, beta []Sample) bool {
	if len(alpha) != len(beta) {
		return false
	}
	for i := range alpha {
		if len(alpha[i].Alleles) != len(beta[i].Alleles) {
			return false
		}
		if len(alpha[i].Phase) != len(beta[i].Phase) {
			return false
		}
		if len(alpha[i].FormatData) != len(beta[i].FormatData) {
			return false
		}

		for j := range alpha[i].Alleles {
			if alpha[i].Alleles[j] != beta[i].Alleles[j] {
				return false
			}
		}

		for j := range alpha[i].Phase {
			if alpha[i].Phase[j] != beta[i].Phase[j] {
				return false
			}
		}

		for j := range alpha[i].FormatData {
			if alpha[i].FormatData[j] != beta[i].FormatData[j] {
				return false
			}
		}
	}
	return true
}

//AllEqual returns true if each Vcf in a slice of Vcf structs contain identical information.
func AllEqual(alpha []Vcf, beta []Vcf) bool {
	if len(alpha) != len(beta) {
		if debugMode > 0 {
			log.Printf("AllEqual is false. len(a): %v. len(b): %v.\n", len(alpha), len(beta))
		}
		return false
	}
	for i := 0; i < len(alpha); i++ {
		if !isEqual(alpha[i], beta[i]) {
			if debugMode > 0 {
				log.Printf("AllEqual is false. Variants at index %v do not match.\n", i)
			}
			return false
		}
	}
	return true
}
