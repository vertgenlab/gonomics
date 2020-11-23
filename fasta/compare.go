package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"sort"
	"strings"
)

func compareName(alpha *Fasta, beta *Fasta) int {
	return strings.Compare(alpha.Name, beta.Name)
}

func compareSeq(alpha *Fasta, beta *Fasta) int {
	return dna.CompareSeqsCaseSensitive(alpha.Seq, beta.Seq)
}

func compareSeqIgnoreCase(alpha *Fasta, beta *Fasta) int {
	return dna.CompareSeqsIgnoreCase(alpha.Seq, beta.Seq)
}

//IsEqual returns true if two input Fasta structs have an equal name and sequence.
func IsEqual(alpha *Fasta, beta *Fasta) bool {
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
		if !IsEqual(alpha[idx], beta[idx]) {
			return false
		}
	}
	return true
}

//AllAreEqual returns true if every entry in a slice of Fasta structs passes IsEqual. Sensitive to order in the slice.
func AllAreEqual(alpha []*Fasta, beta []*Fasta) bool {
	return allEqual(alpha, beta, false)
}

//AllAreEqualIgnoreOrder returns true if every entry in a slice of Fasta structs passes IsEqual. Not sensitive to order in the slice.
func AllAreEqualIgnoreOrder(alpha []*Fasta, beta []*Fasta) bool {
	return allEqual(alpha, beta, true)
}

func SortByName(seqs []*Fasta) {
	sort.Slice(seqs, func(i, j int) bool { return compareName(seqs[i], seqs[j]) == -1 })
}

func SortBySeq(seqs []*Fasta) {
	sort.Slice(seqs, func(i, j int) bool { return compareSeqIgnoreCase(seqs[i], seqs[j]) == -1 })
}

//Note: for QuerySeq, RefPosToAlnPos is probably not required if you are using an assembly fasta as the reference, but if you are querying from alignment Fas, you'll want to get the alnIndex before calling this function

//QuerySeq takes in a slice of fastas and a position (name and index) and returns true if a query sequence of bases matches the fasta at this position.
func QuerySeq(records []*Fasta, chr string, index int, query []dna.Base) bool {
	chrIndex := GetChromIndex(records, chr)
	return dna.CompareSeqsIgnoreCaseAndGaps(query, records[chrIndex].Seq[index:index+len(query)]) != 0
}




