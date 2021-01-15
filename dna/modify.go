package dna

import (
	"errors"
)

var (
	ErrInvalidInterval          = errors.New("deletion interval is not valid")
	ErrInvalidInsertionPosition = errors.New("insertion position is not valid")
	ErrUnrecognizedBase         = errors.New("input base was not recognized")
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

// ComplementSingleBase returns the nucleotide complementary to the input base.
func ComplementSingleBase(b Base) (Base, error) {
	switch b {
	case A:
		return T, nil
	case C:
		return G, nil
	case G:
		return C, nil
	case T:
		return A, nil
	case N:
		return N, nil
	case LowerA:
		return LowerT, nil
	case LowerC:
		return LowerG, nil
	case LowerG:
		return LowerC, nil
	case LowerT:
		return LowerA, nil
	case LowerN:
		return LowerN, nil
	case Gap:
		return Gap, nil
	case Dot:
		return Dot, nil
	case Nil:
		return Nil, nil
	default:
		return b, ErrUnrecognizedBase
	}
}

// ReverseComplement reverses a sequence of bases and complements each base.
// Used to switch strands and maintain 5' -> 3' orientation.
func ReverseComplement(bases []Base) error {
	var err error
	for i, j := 0, len(bases)-1; i <= j; i, j = i+1, j-1 {
		bases[i], bases[j], err = complementSwap(bases[i], bases[j])
		if err != nil {
			break
		}
	}
	return err
}

// Complement all bases in a sequence of bases.
func Complement(bases []Base) error {
	var err error
	for i := range bases {
		bases[i], err = ComplementSingleBase(bases[i])
		if err != nil {
			break
		}
	}
	return err
}

// complementSwap complements the input bases and swaps their positions in the return.
func complementSwap(alpha, beta Base) (Base, Base, error) {
	var err error
	alpha, err = ComplementSingleBase(alpha)
	if err != nil {
		return beta, alpha, err
	}
	beta, err = ComplementSingleBase(beta)
	return beta, alpha, err
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

// Delete removes bases from and sequence of bases.
// all base positions are zero based and left closed, right open.
func Delete(seq []Base, delStart int, delEnd int) ([]Base, error) {
	if delStart >= delEnd || delStart < 0 || delEnd > len(seq) {
		return nil, ErrInvalidInterval
	}
	return append(seq[:delStart], seq[delEnd:]...), nil
}

// Insert adds bases to a sequence of bases.
// base position is zero-based, insertion happens before specified base
// giving the length of the sequence puts the insertion at the end.
func Insert(seq []Base, insPos int, insSeq []Base) ([]Base, error) {
	if insPos < 0 || insPos > len(seq) {
		return nil, ErrInvalidInsertionPosition
	}
	return append(seq[:insPos], append(insSeq, seq[insPos:]...)...), nil
}

// Replace performs both a deletion and an insertion,
// replacing the input interval with the input insSeq.
// all base positions are zero based and left closed, right open.
func Replace(seq []Base, start int, end int, insSeq []Base) []Base {
	return append(seq[:start], append(insSeq, seq[end:]...)...)
}
