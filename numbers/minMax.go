package numbers

import (
	"log"

	"golang.org/x/exp/constraints"
)

// Max returns the maximum of two input values of an ordered type.
func Max[E constraints.Ordered](a, b E) E {
	if a > b {
		return a
	}
	return b
}

// Min returns the minimum of two input values of an ordered type.
func Min[E constraints.Ordered](a, b E) E {
	if a < b {
		return a
	}
	return b
}

// MaxMany returns the maximum of any number of input values of an ordered type.
func MaxMany[E constraints.Ordered](s ...E) E {
	if len(s) == 0 {
		log.Panic("empty argument to max function")
	}
	var max E = s[0]
	for i := 1; i < len(s); i++ {
		if s[i] > max {
			max = s[i]
		}
	}
	return max
}

// MinMany returns the minimum of any number of input values of an ordered type.
func MinMany[E constraints.Ordered](s ...E) E {
	if len(s) == 0 {
		log.Panic("empty argument to min function")
	}
	var min E = s[0]
	for i := 1; i < len(s); i++ {
		if s[i] < min {
			min = s[i]
		}
	}
	return min
}
