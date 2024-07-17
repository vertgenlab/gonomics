package cigar

// AddCigar adds a cigar struct to the end of a slice, but is smart about checking
// to see if the addition is the same operation that is present at the end of the
// existing slice and just increasing the run length of that entry.
func AddCigar(cigs []Cigar, newCig Cigar) []Cigar {
	if len(cigs) == 0 {
		cigs = append(cigs, newCig)
	} else if cigs[len(cigs)-1].Op == newCig.Op {
		cigs[len(cigs)-1].RunLength += newCig.RunLength
	} else {
		cigs = append(cigs, newCig)
	}
	return cigs
}

// CatCigar cats two cigars together, but is smart about checking to see if the
// last element of the first slice is the same operation as the first element
// of the second slice.  In this case it will compress that struct into a single
// element with a longer run length.
func CatCigar(cigs []Cigar, newCigs []Cigar) []Cigar {
	if len(newCigs) == 0 || newCigs == nil {
		return cigs
	} else if len(cigs) == 0 {
		return newCigs
	} else {
		cigs = AddCigar(cigs, newCigs[0])
		cigs = append(cigs, newCigs[1:]...)
		return cigs
	}
}

// ReverseBytesCigar cigar will reverse the order of a cigar slice. Typically performed after matrix traceback
func ReverseBytesCigar(alpha []Cigar) {
	var i, off int
	for i = len(alpha)/2 - 1; i >= 0; i-- {
		off = len(alpha) - 1 - i
		alpha[i], alpha[off] = alpha[off], alpha[i]
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
