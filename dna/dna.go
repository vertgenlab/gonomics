package dna

import (
	"github.com/vertgenlab/gonomics/common"
)

type Base byte

// TODO: change these names so that all variables can be seen
// from outside.  Maybe BaseA, Basea, etc
const (
	A   Base = 0
	C   Base = 1
	G   Base = 2
	T   Base = 3
	N   Base = 4
	a   Base = 5
	c   Base = 6
	g   Base = 7
	t   Base = 8
	n   Base = 9
	Gap Base = 10
)

func ToUpper(b Base) Base {
	if b == Gap {
		return b
	} else if b > N {
		return b - 5
	} else {
		return b
	}
}

func ToLower(b Base) Base {
	if b == Gap {
		return b
	} else if b < a {
		return b + 5
	} else {
		return b
	}
}

// start is closed, end is open, both are zero-based
func RangeToUpper(bases []Base, start int, end int) {
	for i := start; i < end; i++ {
		bases[i] = ToUpper(bases[i])
	}
}

// start is closed, end is open, both are zero-based
func RangeToLower(bases []Base, start int, end int) {
	for i := start; i < end; i++ {
		bases[i] = ToLower(bases[i])
	}
}

func AllToUpper(bases []Base) {
	RangeToUpper(bases, 0, len(bases))
}

func AllToLower(bases []Base) {
	RangeToLower(bases, 0, len(bases))
}

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

func compareSeqs(alpha []Base, beta []Base, ignoreCase bool) int {
	var res int
	stop := common.Min(len(alpha), len(beta))
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
