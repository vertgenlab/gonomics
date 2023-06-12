package fastq

import (
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
)

const debugMode = 0 //set debugMode to 1 to enable prints

// AllAreEqual returns true if each Fastq entry in a slice of Fastq structs contains identical information, false otherwise.
func AllAreEqual(alpha []Fastq, beta []Fastq) bool {
	if len(alpha) != len(beta) {
		if debugMode > 0 {
			log.Printf("AllAreEqual is false. len(alpha): %d. len(beta): %d.\n", len(alpha), len(beta))
		}
		return false
	}
	for i := range alpha {
		if !IsEqual(alpha[i], beta[i]) {
			if debugMode > 0 {
				log.Printf("AllAreEqual is false. Entries at index %d do not match.", i)
			}
			return false
		}
	}
	return true
}

// IsEqual returns true if two input Fastq structs contain identical information, false otherwise.
func IsEqual(alpha Fastq, beta Fastq) bool {
	if CompareName(alpha, beta) != 0 {
		if debugMode == 1 {
			log.Printf("IsEqual: Name is not equal. Alpha: %s. Beta: %s.\n", alpha.Name, beta.Name)
		}
		return false
	}
	if CompareSeq(alpha, beta) != 0 {
		if debugMode == 1 {
			log.Printf("IsEqual: Seq is not equal. Alpha: %s. Beta: %s.\n", dna.BasesToString(alpha.Seq), dna.BasesToString(beta.Seq))
		}
		return false
	}
	if CompareQuals(alpha, beta) != 0 {
		return false
	}
	return true
}

// CompareName compares two Fastq names for sorting or equality testing.
func CompareName(alpha Fastq, beta Fastq) int {
	return strings.Compare(alpha.Name, beta.Name)
}

// CompareSeq compares two Fastq sequences for sorting or equality testing.
func CompareSeq(alpha Fastq, beta Fastq) int {
	return dna.CompareSeqsCaseSensitive(alpha.Seq, beta.Seq)
}

// CompareQuals compares two Fastq quality score slices for sorting or equality testing.
func CompareQuals(alpha Fastq, beta Fastq) int {
	var result int
	stop := min(len(alpha.Qual), len(beta.Qual))
	for i := 0; i < stop; i++ {
		result = CompareQual(alpha.Qual[i], beta.Qual[i])
		if result != 0 {
			return result
		}
	}
	if len(alpha.Qual) < len(beta.Qual) {
		return -1
	} else if len(alpha.Qual) > len(beta.Qual) {
		return 1
	}
	return 0
}

// CompareQual compares two Fastq quality scores for sorting or equality testing.
func CompareQual(alpha uint8, beta uint8) int {
	if alpha < beta {
		return -1
	} else if alpha > beta {
		return 1
	}
	return 0
}

// min helper function to avoid dependencies.
func min(a int, b int) int {
	if a <= b {
		return a
	} else {
		return b
	}
}

// PairedEndIsEqual returns true if two input PairedEnd structs contain identical information, false otherwise.
func PairedEndIsEqual(alpha PairedEnd, beta PairedEnd) bool {
	if !IsEqual(alpha.Fwd, beta.Fwd) {
		if debugMode == 1 {
			log.Printf("PairedEndIsEqual: Forward is not equal.\n")
		}
		return false
	}
	if !IsEqual(alpha.Rev, beta.Rev) {
		if debugMode == 1 {
			log.Printf("PairedEndIsEqual: Reverse is not equal.\n")
		}
		return false
	}
	return true
}

// SingleCellIsEqual returns true if two input SingleCellPair structs are equal, false otherwise.
func SingleCellIsEqual(alpha SingleCellPair, beta SingleCellPair) bool {
	if !PairedEndIsEqual(alpha.Reads, beta.Reads) {
		if debugMode == 1 {
			log.Printf("SingleCellIsEqual: PairedEnds were not equal.\n")
		}
		return false
	}
	if dna.CompareSeqsCaseSensitive(alpha.Bx, beta.Bx) != 0 {
		if debugMode == 1 {
			log.Printf("SingleCellIsEqual: Barcodes were not equal. Alpha: %s. Beta: %s.\n", dna.BasesToString(alpha.Bx), dna.BasesToString(beta.Bx))
		}
		return false
	}
	if dna.CompareSeqsCaseSensitive(alpha.Umi, beta.Umi) != 0 {
		if debugMode == 1 {
			log.Printf("SingleCellIsEqual: Umis were not equal. Alpha: %s. Beta: %s.\n", dna.BasesToString(alpha.Umi), dna.BasesToString(beta.Umi))
		}
		return false
	}
	return true
}
