package cigar

// Append merges a cigar into a slice, checking for consecutive cigars with the same operation.
func Append(alpha []Cigar, beta Cigar) []Cigar {
	if len(alpha) > 0 && alpha[len(alpha)-1].Op == beta.Op {
		alpha[len(alpha)-1].RunLength += beta.RunLength
	} else {
		alpha = append(alpha, beta)
	}
	return alpha
}

// Concat combines two cigar slices, checking consecutive elements with the same operation.
func Concat(alpha []Cigar, beta []Cigar) []Cigar {
	if len(alpha) == 0 {
		return beta
	}
	if len(beta) > 0 {
		alpha = Append(alpha, beta[0])
		beta = beta[1:]
	}
	return append(alpha, beta...)
}

// AppendSoftClips adds soft clips to the beginning and/or end of a CIGAR string to match a given read length.
func AppendSoftClips(front, lengthOfRead int, cigars []Cigar) []Cigar {
	var currRunLen = QueryLength(cigars) // Fixed variable name
	if front == 0 && currRunLen >= lengthOfRead {
		return cigars
	}
	// Pre-allocate slice for efficiency (estimate maximum size)
	var answer []Cigar = make([]Cigar, 0, len(cigars)+2)
	if front > 0 {
		answer = append(answer, Cigar{RunLength: front, Op: SoftClip})
	}
	if front+currRunLen < lengthOfRead {
		answer = append(append(answer, cigars...), Cigar{RunLength: lengthOfRead - front - currRunLen, Op: SoftClip})
	}
	return answer
}

// ReverseCigar reverses the order of a cigar slice in place.
func ReverseCigar(cigars []Cigar) {
	size := len(cigars)
	for i, j := 0, size-1; i < size/2; i, j = i+1, j-1 {
		cigars[i], cigars[j] = cigars[j], cigars[i]
	}
}

// ToUint32 encodes cigar op and runlen as a uint32 defined by runlength(op)<<4|op.
func ToUint32(c Cigar) uint32 {
	return uint32(c.RunLength)<<4 | Uint32Table[c.Op] // move 4 bits to the left followed by bitwise OR with op.
}

// TODO: Move TripleMaxTrace() and TripleMaxTraceExtended() to align package and replace
// TripleMaxTrace will trace smith-waterman matrix alignment and return one of 3 cigar Op's.
// M: matches or mismatches, I: insertions, D: for deletions.
func TripleMaxTrace(a int64, b int64, c int64) (int64, byte) {
	if a >= b && a >= c {
		return a, Match
	} else if b >= c {
		return b, Insertion
	} else {
		return c, Deletion
	}
}

// TripleMaxTraceExtended is an expanded version of TripleMaxTrace which trace will return either '=' or 'X' where 'M'(s) are found.
func TripleMaxTraceExtended(prev int64, a int64, b int64, c int64) (int64, byte) {
	if a >= b && a >= c {
		if a > prev {
			return a, Equal
		}
		return a, Mismatch
	}
	if b >= c {
		return b, Insertion
	}
	return c, Deletion
}

// func traceBack(prev int64, a int64, b int64, c int64) (int64, byte) {
//     if a >= b && a >= c {
//         if a > prev {
//             return a, cigar.Equal
//         }
//         return a, cigar.Mismatch
//     }
//     if b >= c {
//         return b, cigar.Insertion
//     }
//     return c, cigar.Deletion
// }
