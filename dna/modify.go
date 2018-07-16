package dna

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
)

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
