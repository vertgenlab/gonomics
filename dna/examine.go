package dna

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
