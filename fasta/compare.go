package fasta

import (
	"sort"

	"github.com/vertgenlab/gonomics/dna"
)

// IsEqual returns true if two input Fasta structs have an equal name and sequence.
func IsEqual(alpha Fasta, beta Fasta) bool {
	if alpha.Name == beta.Name && dna.CompareSeqsCaseSensitive(alpha.Seq, beta.Seq) == 0 {
		return true
	} else {
		return false
	}
}

// allEqual determines if two slices of fasta records are equivalent.
func allEqual(alpha []Fasta, beta []Fasta, ignoreOrder bool) bool {
	if len(alpha) != len(beta) {
		return false
	}
	if ignoreOrder {
		SortByName(alpha)
		SortByName(beta)
	}
	for i := range alpha {
		if !IsEqual(alpha[i], beta[i]) {
			return false
		}
	}
	return true
}

// AllAreEqual returns true if every entry in a slice of Fasta structs passes IsEqual. Sensitive to order in the slice.
func AllAreEqual(alpha []Fasta, beta []Fasta) bool {
	return allEqual(alpha, beta, false)
}

// AllAreEqualIgnoreOrder returns true if every entry in a slice of Fasta structs passes IsEqual. Not sensitive to order in the slice.
func AllAreEqualIgnoreOrder(alpha []Fasta, beta []Fasta) bool {
	return allEqual(alpha, beta, true)
}

// SortByName sorts fasta records lexicographically.
func SortByName(seqs []Fasta) {
	sort.Slice(seqs, func(i, j int) bool { return seqs[i].Name < seqs[j].Name })
}

// SortBySeq sorts fasta records by sequence.
func SortBySeq(seqs []Fasta) {
	sort.Slice(seqs, func(i, j int) bool { return dna.CompareSeqsIgnoreCase(seqs[i].Seq, seqs[j].Seq) == -1 })
}
