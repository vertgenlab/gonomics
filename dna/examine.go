package dna

import (
	"fmt"
	"log"
	"github.com/vertgenlab/gonomics/common"
)

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
		case a:
			aCount++
		case c:
			cCount++
		case g:
			gCount++
		case t:
			tCount++
		case n:
			nCount++
		case Gap:
			gapCount++
		}
	}
	return ACount, CCount, GCount, TCount, NCount, aCount, cCount, gCount, tCount, NCount, gapCount
}

func CountMask(seq []Base) (unmaskedCount int, maskedCount int, gapCount int) {
	ACount, CCount, GCount, TCount, NCount, aCount, cCount, gCount, tCount, nCount, gapCount := Count(seq)
	unmaskedCount = ACount + CCount + GCount + TCount + NCount
	maskedCount = aCount + cCount + gCount + tCount + nCount
	return unmaskedCount, maskedCount, gapCount
}

func CountGaps(seq []Base) int {
	_, _, gapCount := CountMask(seq)
	return gapCount
}

func BaseDist(a Base, b Base) int {
	if a == b {
		return 0
	}
	return 1
}

func Dist(a []Base, b []Base) int {
	if len(a) != len(b) {
		log.Fatalf("Seqs must have the same length to calculate distance.\n")
	}
	var sum int = 0
	for i := 0; i < len(a); i++ {
		sum = sum + BaseDist(a[i], b[i])
	}
	return sum
}

func IsLower(b Base) bool {
	switch b{
	case a:
		return true
	case g:
		return true
	case c:
		return true
	case t:
		return true
	default:
		return false
	}
}

func DefineBase(b Base) bool {
	switch b {
	case A:
		return true
	case C:
		return true
	case G:
		return true
	case T:
		return true
	case N:
		return false
	case a:
		return true
	case c:
		return true
	case g:
		return true
	case t:
		return true
	case n:
		return false
	case Gap:
		return false
	default:
		common.ExitIfError(fmt.Errorf("Error: trying to examine unknown base %d", b))
		return false
	}
}
