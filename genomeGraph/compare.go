package genomeGraph

import (
	"sort"

	"github.com/vertgenlab/gonomics/dna"
)

func compareId(alpha *Node, beta *Node) int {
	if alpha.Id < beta.Id {
		return 1
	} else if beta.Id < alpha.Id {
		return -1
	} else {
		return 0
	}
}

func compareSeq(alpha *Node, beta *Node) int {
	return dna.CompareSeqsCaseSensitive(alpha.Seq, beta.Seq)
}

func compareSeqIgnoreCase(alpha *Node, beta *Node) int {
	return dna.CompareSeqsIgnoreCase(alpha.Seq, beta.Seq)
}

func isEqual(alpha *Node, beta *Node) bool {
	if compareId(alpha, beta) == 0 && compareSeq(alpha, beta) == 0 {
		return true
	} else {
		return false
	}
}

func allEqual(alpha []*Node, beta []*Node, ignoreOrder bool) bool {
	if len(alpha) != len(beta) {
		return false
	}
	if ignoreOrder {
		SortById(alpha)
		SortById(beta)
	}
	for idx := range alpha {
		if !isEqual(alpha[idx], beta[idx]) {
			return false
		}
	}
	return true
}

func AllAreEqual(alpha []*Node, beta []*Node) bool {
	return allEqual(alpha, beta, false)
}

func AllAreEqualIgnoreOrder(alpha []*Node, beta []*Node) bool {
	return allEqual(alpha, beta, true)
}

func SortById(seqs []*Node) {
	sort.Slice(seqs, func(i, j int) bool { return compareId(seqs[i], seqs[j]) == -1 })
}

func SortBySeq(seqs []*Node) {
	sort.Slice(seqs, func(i, j int) bool { return compareSeqIgnoreCase(seqs[i], seqs[j]) == -1 })
}
