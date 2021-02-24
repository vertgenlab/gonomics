package fasta

import (
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
