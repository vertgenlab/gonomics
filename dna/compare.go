package dna

import (
	"github.com/vertgenlab/gonomics/numbers"
	"log"
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

func compareSeqsIgnoreGaps(alpha []Base, beta []Base, ignoreCase bool) int {
	var res, i, j int

	for i, j = 0, 0; i < len(alpha) && j < len(beta); i, j = i+1, j+1 {
		for ; i < len(alpha) && alpha[i] == Gap; i++ {
		}
		for ; j < len(beta) && beta[j] == Gap; j++ {
		}
		if i < len(alpha) && j < len(beta) {
			res = compareBases(alpha[i], beta[j], ignoreCase)
			if res != 0 {
				return res
			}
		}
	}
	for ; i < len(alpha) && alpha[i] == Gap; i++ {
	}
	for ; j < len(beta) && beta[j] == Gap; j++ {
	}
	if i == len(alpha) && j == len(beta) {
		return 0
	} else if i == len(alpha) {
		return -1
	} else if j == len(beta) {
		return 1
	} else {
		log.Fatal("Error: unexpectedly got to else in if statement\n")
		return 0
	}
}

func compareSeqs(alpha []Base, beta []Base, ignoreCase bool) int {
	var res int
	stop := numbers.Min(len(alpha), len(beta))
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

func CompareSeqsIgnoreCaseAndGaps(alpha []Base, beta []Base) int {
	return compareSeqsIgnoreGaps(alpha, beta, true)
}

func CompareSeqsCaseSensitiveIgnoreGaps(alpha []Base, beta []Base) int {
	return compareSeqsIgnoreGaps(alpha, beta, false)
}

func compare2DSeqs(alpha [][]Base, beta [][]Base, ignoreCase bool, ignoreGaps bool) int {
	var res int
	stop := numbers.Min(len(alpha), len(beta))
	for i := 0; i < stop; i++ {
		if ignoreGaps {
			res = compareSeqsIgnoreGaps(alpha[i], beta[i], ignoreCase)
		} else {
			res = compareSeqs(alpha[i], beta[i], ignoreCase)
		}
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

func CompareTwoDSeqsIgnoreCase(alpha [][]Base, beta [][]Base) int {
	return compare2DSeqs(alpha, beta, true, false)
}

func CompareTwoDSeqsCaseSensitive(alpha [][]Base, beta [][]Base) int {
	return compare2DSeqs(alpha, beta, false, false)
}

func CompareTwoDSeqsIgnoreCaseAndGaps(alpha [][]Base, beta [][]Base) int {
	return compare2DSeqs(alpha, beta, true, true)
}

func CompareTwoDSeqsCaseSensitiveIgnoreGaps(alpha [][]Base, beta [][]Base) int {
	return compare2DSeqs(alpha, beta, false, true)
}
