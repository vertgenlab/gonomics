package fasta

import (
	"github.com/vertgenlab/gonomics/chromInfo"
	"log"
	"strings"
)

// IsFasta returns true if the input filename has a fasta file extension.
// Input filename may have a .gz suffix.
func IsFasta(filename string) bool {
	if strings.HasSuffix(filename, ".fasta") ||
		strings.HasSuffix(filename, ".fa") ||
		strings.HasSuffix(filename, ".fasta.gz") ||
		strings.HasSuffix(filename, ".fa.gz") {
		return true
	} else {
		return false
	}
}

// ToChromInfo converts a []Fasta into a []ChromInfo. Useful for applications
// that do not require the entire fasta sequence to be kept in memory, but just
// the name, size, and order of fasta records.
func ToChromInfo(records []Fasta) []chromInfo.ChromInfo {
	answer := make([]chromInfo.ChromInfo, len(records))
	for i := range records {
		answer[i] = chromInfo.ChromInfo{Name: records[i].Name, Size: len(records[i].Seq), Order: i}
	}
	return answer
}

// GetChromIndex returns the index of a Fasta record whose name matches an input string
func GetChromIndex(f []Fasta, chrom string) int {
	for i := range f {
		if f[i].Name == chrom {
			return i
		}
	}
	log.Fatalf("Error. Chromosome name: %s not found in fasta.", chrom)
	return -1
}
