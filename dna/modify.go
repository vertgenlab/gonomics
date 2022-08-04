package dna

import (
	"log"
)

// ToUpper changes the input base to uppercase.
func ToUpper(b Base) Base {
	switch b {
	case LowerA:
		return A
	case LowerC:
		return C
	case LowerG:
		return G
	case LowerT:
		return T
	case LowerN:
		return N
	default:
		return b
	}
}

// ToLower changes the input base to lowercase.
func ToLower(b Base) Base {
	switch b {
	case A:
		return LowerA
	case C:
		return LowerC
	case G:
		return LowerG
	case T:
		return LowerT
	case N:
		return LowerN
	default:
		return b
	}
}

// RangeToUpper changes the bases in a set range to uppercase.
// start is closed, end is open, both are zero-based.
func RangeToUpper(bases []Base, start int, end int) {
	for i := start; i < end; i++ {
		bases[i] = ToUpper(bases[i])
	}
}

// RangeToLower changes the bases in a set range to lowercase.
// start is closed, end is open, both are zero-based.
func RangeToLower(bases []Base, start int, end int) {
	for i := start; i < end; i++ {
		bases[i] = ToLower(bases[i])
	}
}

// AllToUpper changes all bases in a sequence to uppercase.
func AllToUpper(bases []Base) {
	RangeToUpper(bases, 0, len(bases))
}

// AllToLower changes all bases in a sequence to lowercase.
func AllToLower(bases []Base) {
	RangeToLower(bases, 0, len(bases))
}

// complementArray is an efficient lookup for the complement of the input base.
// intended to remain as a private array to help the Complement functions.
// panics if value input is not a valid Base.
var complementArray = []Base{T, G, C, A, N, LowerT, LowerG, LowerC, LowerA, LowerN, Gap, Dot, Nil}

// ComplementSingleBase returns the nucleotide complementary to the input base.
func ComplementSingleBase(b Base) Base {
	switch b {
	case A:
		return T
	case C:
		return G
	case G:
		return C
	case T:
		return A
	case N:
		return N
	case LowerA:
		return LowerT
	case LowerC:
		return LowerG
	case LowerG:
		return LowerC
	case LowerT:
		return LowerA
	case LowerN:
		return LowerN
	case Gap:
		return Gap
	case Dot:
		return Dot
	case Nil:
		return Nil
	default:
		log.Panicf("unrecognized base %v", b)
		return Nil
	}
}

//Reverse an input sequence of bases. Does not maintain 5' -> 3' orientation.
func Reverse(bases []Base) {
	for i, j := 0, len(bases)-1; i <= j; i, j = i+1, j-1 {
		bases[i], bases[j] = bases[j], bases[i]
	}
}

// ReverseComplement reverses a sequence of bases and complements each base.
// Used to switch strands and maintain 5' -> 3' orientation.
func ReverseComplement(bases []Base) {
	for i, j := 0, len(bases)-1; i <= j; i, j = i+1, j-1 {
		bases[i], bases[j] = complementArray[bases[j]], complementArray[bases[i]]
	}
}

// Complement all bases in a sequence of bases.
func Complement(bases []Base) {
	for i := range bases {
		bases[i] = complementArray[bases[i]]
	}
}

// RemoveGaps returns a sequence of bases with no gaps.
func RemoveGaps(bases []Base) []Base {
	return RemoveBase(bases, Gap)
}

// RemoveBase returns a sequence of bases without any of the designated base.
func RemoveBase(bases []Base, baseToRemove Base) []Base {
	ans := make([]Base, 0, len(bases))
	for i := range bases {
		if bases[i] != baseToRemove {
			ans = append(ans, bases[i])
		}
	}
	return ans
}

// Delete removes bases from a sequence of bases.
// all base positions are zero based and left closed, right open.
func Delete(seq []Base, delStart int, delEnd int) []Base {
	if delStart >= delEnd || delStart < 0 || delEnd > len(seq) {
		log.Panicf("a deletion on sequence of length %d with start of %d and length of %d is invalid", len(seq), delStart, delEnd)
	}
	return append(seq[:delStart], seq[delEnd:]...)
}

// Insert adds bases to a sequence of bases.
// base position is zero-based, insertion happens before specified base
// giving the length of the sequence puts the insertion at the end.
func Insert(seq []Base, insPos int, insSeq []Base) []Base {
	if insPos < 0 || insPos > len(seq) {
		log.Panicf("an insertion on sequence of length %d with start of %d is invalid", len(seq), insPos)
	}
	return append(seq[:insPos], append(insSeq, seq[insPos:]...)...)
}

// Replace performs both a deletion and an insertion,
// replacing the input interval with the input insSeq.
// all base positions are zero based and left closed, right open.
func Replace(seq []Base, start int, end int, insSeq []Base) []Base {
	return append(seq[:start], append(insSeq, seq[end:]...)...)
}
