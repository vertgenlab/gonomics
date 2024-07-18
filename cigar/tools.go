package cigar

// AddCigar adds a cigar struct to the end of a slice, but is smart about checking
// to see if the addition is the same operation that is present at the end of the
// existing slice and just increasing the run length of that entry.
func AddCigar(alpha []Cigar, beta Cigar) []Cigar {
	if len(alpha) > 0 && alpha[len(alpha)-1].Op == beta.Op {
		alpha[len(alpha)-1].RunLength += beta.RunLength
	} else {
		alpha = append(alpha, beta)
	}
	return alpha
}

// CatCigar cats two cigars together, but is smart about checking to see if the
// last element of the first slice is the same operation as the first element
// of the second slice.  In this case it will compress that struct into a single
// element with a longer run length.
func CatCigar(alpha []Cigar, beta []Cigar) []Cigar {
	if len(beta) == 0 {
		return alpha
	}

	size := len(alpha)
	if size == 0 {
		return beta
	}
	// Merge the first new cigar if possible
	if alpha[size-1].Op == beta[0].Op {
		alpha[size-1].RunLength += beta[0].RunLength
		beta = beta[1:] // Remove the merged cigar from beta
	}
	// Append the remaining new cigars
	return append(alpha, beta...)
}

// ReverseCigar reverses the order of a cigar slice in place.
func ReverseCigar(cigars []Cigar) {
	for i, j := 0, len(cigars)-1; i < j; i, j = i+1, j-1 {
		cigars[i], cigars[j] = cigars[j], cigars[i]
	}
}

// MatrixTrace will trace smith-waterman matrix alignment and return one of 3 cigar Op's.
// M: matches or mismatches, I: insertions, D: for deletions.
func MatrixTrace(a int64, b int64, c int64) (int64, byte) {
	if a >= b && a >= c {
		return a, Match
	} else if b >= c {
		return b, Insertion
	} else {
		return c, Deletion
	}
}

// TripleMaxTrace is an expanded version of MatrixTrace which will return either '=' or 'X' where 'M'(s) are found.
func TripleMaxTrace(prev int64, a int64, b int64, c int64) (int64, byte) {
	if a >= b && a >= c {
		if a > prev {
			return a, Equal
		} else {
			return a, Mismatch
		}
	} else if b >= c {
		return b, Insertion
	} else {
		return c, Deletion
	}
}
