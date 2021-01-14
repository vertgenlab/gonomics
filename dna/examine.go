package dna

import (
	"errors"
)

var (
	ErrInputSeqsDiffLen = errors.New("input sequences are not the same length")
)

// Count returns the number of each base present in the input sequence.
func Count(seq []Base) (ACount int, CCount int, GCount int, TCount int, NCount int, aCount int, cCount int, gCount int, tCount int, nCount int, gapCount int) {
	ACount, CCount, GCount, TCount, NCount, aCount, cCount, gCount, tCount, nCount, gapCount = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	for _, b := range seq {
		switch b {
		case A:
			ACount++
		case C:
			CCount++
		case G:
			GCount++
		case T:
			TCount++
		case N:
			NCount++
		case LowerA:
			aCount++
		case LowerC:
			cCount++
		case LowerG:
			gCount++
		case LowerT:
			tCount++
		case LowerN:
			nCount++
		case Gap:
			gapCount++
		}
	}
	return ACount, CCount, GCount, TCount, NCount, aCount, cCount, gCount, tCount, NCount, gapCount
}

// CountMask returns the number of bases that are masked/unmasked (lowercase/uppercase) in the input sequence.
func CountMask(seq []Base) (unmaskedCount int, maskedCount int, gapCount int) {
	ACount, CCount, GCount, TCount, NCount, aCount, cCount, gCount, tCount, nCount, gapCount := Count(seq)
	unmaskedCount = ACount + CCount + GCount + TCount + NCount
	maskedCount = aCount + cCount + gCount + tCount + nCount
	return unmaskedCount, maskedCount, gapCount
}

// CountGaps returns the number of gaps present in the input sequence.
func CountGaps(seq []Base) int {
	var gapCount int
	for _, val := range seq {
		if val == Gap {
			gapCount++
		}
	}
	return gapCount
}

// baseDist is a helper function for Dist that returns 1 if input bases do not match.
func baseDist(a Base, b Base) int {
	if a == b {
		return 0
	}
	return 1
}

// Dist returns the number of bases that do not match between the input sequences.
// Input sequences must be the sample length.
func Dist(a []Base, b []Base) (int, error) {
	if len(a) != len(b) {
		return 0, ErrInputSeqsDiffLen
	}
	var sum int
	for i := range a {
		sum = sum + baseDist(a[i], b[i])
	}
	return sum, nil
}

// IsLower returns true if the input base is lowercase.
func IsLower(b Base) bool {
	switch b {
	case LowerA:
		return true
	case LowerG:
		return true
	case LowerC:
		return true
	case LowerT:
		return true
	default:
		return false
	}
}

// DefineBase returns false if the input base is an N, Gap, Dot, or Nil.
func DefineBase(b Base) (bool, error) {
	switch b {
	case A:
		return true, nil
	case C:
		return true, nil
	case G:
		return true, nil
	case T:
		return true, nil
	case N:
		return false, nil
	case LowerA:
		return true, nil
	case LowerC:
		return true, nil
	case LowerG:
		return true, nil
	case LowerT:
		return true, nil
	case LowerN:
		return false, nil
	case Gap:
		return false, nil
	case Dot:
		return false, nil
	case Nil:
		return false, nil
	default:
		return false, ErrUnrecognizedBase
	}
}

// CountBase returns the number of the designated base present in the input sequence.
func CountBase(seq []Base, b Base) (int, error) {
	return CountBaseInterval(seq, b, 0, len(seq))
}

// CountBaseInterval returns the number of the designated base present in the input range of the sequence.
func CountBaseInterval(seq []Base, b Base, start int, end int) (int, error) {
	var answer int
	if start < 0 || end > len(seq) {
		return answer, ErrInvalidInterval
	}
	for i := start; i < end; i++ {
		if seq[i] == b {
			answer++
		}
	}
	return answer, nil
}
