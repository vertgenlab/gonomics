package dna

// compareBases returns an integer related to the lexographical order of nucleotides.
// i.e. A < C < a < c < Dot < Gap.
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

// compareSeqsIgnoreCase returns an integer defining the relationship between two input sequences.
// Ignores gaps in the input sequences.
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

	switch {
	case i == len(alpha) && j == len(beta):
		return 0
	case i == len(alpha):
		return -1
	case j == len(beta):
		return 1
	default: // unreachable
		return -1
	}
}

// compareSeqs returns an integer defining the relationship between two input sequences.
func compareSeqs(alpha []Base, beta []Base, ignoreCase bool) int {
	var res int
	stop := min(len(alpha), len(beta))
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

// CompareSeqsIgnoreCase returns an integer defining the relationship between two input sequences.
// 1 if alpha > beta, -1 if beta > alpha, 0 if the sequences are equal.
// Case insensitive.
func CompareSeqsIgnoreCase(alpha []Base, beta []Base) int {
	return compareSeqs(alpha, beta, true)
}

// CompareSeqsCaseSensitive returns an integer defining the relationship between two input sequences.
// 1 if alpha > beta, -1 if beta > alpha, 0 if the sequences are equal.
// Case sensitive.
func CompareSeqsCaseSensitive(alpha []Base, beta []Base) int {
	return compareSeqs(alpha, beta, false)
}

// CompareSeqsIgnoreCaseAndGaps returns an integer defining the relationship between two input sequences.
// 1 if alpha > beta, -1 if beta > alpha, 0 if the sequences are equal.
// Case insensitive. Ignores gaps.
func CompareSeqsIgnoreCaseAndGaps(alpha []Base, beta []Base) int {
	return compareSeqsIgnoreGaps(alpha, beta, true)
}

// CompareSeqsCaseSensitiveIgnoreGaps returns an integer defining the relationship between two input sequences.
// 1 if alpha > beta, -1 if beta > alpha, 0 if the sequences are equal.
// Case sensitive. Ignores gaps.
func CompareSeqsCaseSensitiveIgnoreGaps(alpha []Base, beta []Base) int {
	return compareSeqsIgnoreGaps(alpha, beta, false)
}

// compare2DSeqs returns an integer defining the relationship between two input lists of sequences.
func compare2DSeqs(alpha [][]Base, beta [][]Base, ignoreCase bool, ignoreGaps bool) int {
	var res int
	stop := min(len(alpha), len(beta))
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

// CompareTwoDSeqsIgnoreCase returns an integer defining the relationship between two input lists of sequences.
// 1 if alpha > beta, -1 if beta > alpha, 0 if the sequences are equal.
// Case insensitive.
func CompareTwoDSeqsIgnoreCase(alpha [][]Base, beta [][]Base) int {
	return compare2DSeqs(alpha, beta, true, false)
}

// CompareTwoDSeqsCaseSensitive returns an integer defining the relationship between two input lists of sequences.
// 1 if alpha > beta, -1 if beta > alpha, 0 if the sequences are equal.
// Case sensitive.
func CompareTwoDSeqsCaseSensitive(alpha [][]Base, beta [][]Base) int {
	return compare2DSeqs(alpha, beta, false, false)
}

// CompareTwoDSeqsIgnoreCaseAndGaps returns an integer defining the relationship between two input lists of sequences.
// 1 if alpha > beta, -1 if beta > alpha, 0 if the sequences are equal.
// Case insensitive. Ignores gaps.
func CompareTwoDSeqsIgnoreCaseAndGaps(alpha [][]Base, beta [][]Base) int {
	return compare2DSeqs(alpha, beta, true, true)
}

// CompareTwoDSeqsCaseSensitiveIgnoreGaps returns an integer defining the relationship between two input lists of sequences.
// 1 if alpha > beta, -1 if beta > alpha, 0 if the sequences are equal.
// Case sensitive. Ignores gaps.
func CompareTwoDSeqsCaseSensitiveIgnoreGaps(alpha [][]Base, beta [][]Base) int {
	return compare2DSeqs(alpha, beta, false, true)
}

// IsSeqOfACGT returns true if the input sequences contains only uppercase A/C/G/T.
func IsSeqOfACGT(seq []Base) bool {
	for _, val := range seq {
		if val != A && val != C && val != G && val != T {
			return false
		}
	}
	return true
}

func min(a int, b int) int {
	if a <= b {
		return a
	} else {
		return b
	}
}
