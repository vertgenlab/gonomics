package numbers

import ()

// EqualSliceInt is a function what will check two slices of int if the values are the same.
// Ordering of the slices matter and will only return true if ordering of slice values are the same as well.
func EqualSliceInt(x []int, y []int) bool {
	if len(x) != len(y) {
		return false
	}
	for i := range x {
		if x[i] != y[i] {
			return false
		}
	}
	return true
}
