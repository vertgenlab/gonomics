package numbers

import (
	"golang.org/x/exp/constraints"
	"log"
)

// Max returns the maximum of the input values of an ordered type.
func Max[E constraints.Ordered](s ...E) E {
	if len(s) == 0 {
		log.Panic("empty argument to max function")
	}
	var max E = s[0]
	for i := range s {
		if s[i] > max {
			max = s[i]
		}
	}
	return max
}

//Min returns the minimum of the input values of an ordered type.
func Min[E constraints.Ordered](s ...E) E {
	if len(s) == 0 {
		log.Panic("empty argument to min function")
	}
	var min E = s[0]
	for i := range s {
		if s[i] < min {
			min = s[i]
		}
	}
	return min
}
