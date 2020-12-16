package numbers

import ()

// EqualSliceInt is a function what will check two slices of int if the values are the same.
func EqualSliceInt(x []int, y []int) bool {
	if len(x) != len(y) {
		return false
	}
	for i := 0; i < len(x); i++ {
		if x[i] != y[i] {
			return false
		}
	}
	return true
}
