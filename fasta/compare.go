package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"sort"
	"strings"
)

func compareName(alpha *Fasta, beta *Fasta) int {
	return strings.Compare(alpha.Name, beta.Name)
}

//largest to smallest
func compareLength(alpha *Fasta, beta *Fasta) int {
	if len(alpha.Seq) > len(beta.Seq) {
		return -1
	} else if len(alpha.Seq) < len(beta.Seq) {
		return 1
	} else {
		return 0
	}
}

func compareSeq(alpha *Fasta, beta *Fasta) int {
	return dna.CompareSeqsCaseSensitive(alpha.Seq, beta.Seq)
}

func compareSeqIgnoreCase(alpha *Fasta, beta *Fasta) int {
	return dna.CompareSeqsIgnoreCase(alpha.Seq, beta.Seq)
}

func isEqual(alpha *Fasta, beta *Fasta) bool {
	if compareName(alpha, beta) == 0 && compareSeq(alpha, beta) == 0 {
		return true
	} else {
		return false
	}
}

func allEqual(alpha []*Fasta, beta []*Fasta, ignoreOrder bool) bool {
	if len(alpha) != len(beta) {
		return false
	}
	if ignoreOrder {
		SortByName(alpha)
		SortByName(beta)
	}
	for idx, _ := range alpha {
		if !isEqual(alpha[idx], beta[idx]) {
			return false
		}
	}
	return true
}

func AllAreEqual(alpha []*Fasta, beta []*Fasta) bool {
	return allEqual(alpha, beta, false)
}

func AllAreEqualIgnoreOrder(alpha []*Fasta, beta []*Fasta) bool {
	return allEqual(alpha, beta, true)
}

func SortByName(seqs []*Fasta) {
	sort.Slice(seqs, func(i, j int) bool { return compareName(seqs[i], seqs[j]) == -1 })
}

func SortBySeq(seqs []*Fasta) {
	sort.Slice(seqs, func(i, j int) bool { return compareSeqIgnoreCase(seqs[i], seqs[j]) == -1 })
}

func SortBySeqLen(seqs []*Fasta) {
	sort.Slice(seqs, func(i, j int) bool { return compareLength(seqs[i], seqs[j]) == -1 })
}
