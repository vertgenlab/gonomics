package dna

import (
	"github.com/vertgenlab/gonomics/common"
)

func compareBases(alpha Base, beta Base, ignoreCase bool) int {
	if ignoreCase {
		alpha = ToUpper(alpha)
		beta = ToUpper(beta)
	}
	if alpha < beta {
		return -1
	} else if alpha > beta {
		return 1
	} else {
		return 0
	}
}

func compareSeqs(alpha []Base, beta []Base, ignoreCase bool) int {
	var res int
	stop := common.Min(len(alpha), len(beta))
	for i := 0; i < stop; i++ {
		res = compareBases(alpha[i], beta[i], ignoreCase)
		if res != 0 {
			return res
		}
	}
	if len(alpha) < len(beta) {
		return -1
	} else if len(alpha) > len(beta) {
		return 1
	} else {
		return 0
	}
}

func CompareSeqsIgnoreCase(alpha []Base, beta []Base) int {
	return compareSeqs(alpha, beta, true)
}

func CompareSeqsCaseSensitive(alpha []Base, beta []Base) int {
	return compareSeqs(alpha, beta, false)
}
