package dna

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"log"
)

func ToUpper(b Base) Base {
	if b == Gap || b == Dot {
		return b
	} else if b > N {
		return b - 5
	} else {
		return b
	}
}

func ToLower(b Base) Base {
	if b == Gap || b == Dot {
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

func complement(b Base) Base {
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
	case a:
		return t
	case c:
		return g
	case g:
		return c
	case t:
		return a
	case n:
		return n
	case Gap:
		return Gap
	case Dot:
		return Dot
	default:
		common.ExitIfError(fmt.Errorf("Error: trying to reverse complement an unexpected base %d", b))
		return N
	}
}

func swap(alpha Base, beta Base) (Base, Base) {
	return beta, alpha
}

func ReverseComplement(bases []Base) {
	for i, j := 0, len(bases)-1; i <= j; i, j = i+1, j-1 {
		bases[i], bases[j] = swap(complement(bases[i]), complement(bases[j]))
	}
}

func Complement(bases []Base) {
	for i := 0; i < len(bases); i++ {
		bases[i] = complement(bases[i])
	}
}

func RemoveGaps(bases []Base) []Base {
	var ans []Base
	for i := 0; i < len(bases); i++ {
		if bases[i] != Gap {
			ans = append(ans, bases[i])
		}
	}
	return ans
}

func RemoveBase(bases []Base, b Base) []Base {
	var ans []Base
	for i := 0; i < len(bases); i++ {
		if bases[i] != b {
			ans = append(ans, bases[i])
		}
	}
	return ans
}

// all base positions are zero based and left closed, right open
func Delete(seq []Base, delStart int64, delEnd int64) []Base {
	if delStart >= delEnd || delStart < 0 || delEnd > int64(len(seq)) {
		log.Fatalf("Error: deletion interval from %d to %d is not valid.\n", delStart, delEnd)
	}
	return append(seq[:delStart], seq[delEnd:]...)
}

// base position is zero-based, insertion happens before specified base
// giving the length of the sequence puts the insertion at the end
func Insert(seq []Base, insPos int64, insSeq []Base) []Base {
	if insPos < 0 || insPos > int64(len(seq)) {
		log.Fatalf("Error: insertion location of %d is not valid.\n", insPos)
	}
	return append(seq[:insPos], append(insSeq, seq[insPos:]...)...)
}

// all base positions are zero based and left closed, right open
func Replace(seq []Base, start int64, end int64, insSeq []Base) []Base {
	return append(seq[:start], append(insSeq, seq[end:]...)...)
}
