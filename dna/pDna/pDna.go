package pDna

import "math"

// Float32Base encodes a DNA base as a probability vector using float32 precision.
// Note that gap probabilities can be stored implicitly as 1 - (A + C + G + T).
type Float32Base struct {
	A float32
	C float32
	G float32
	T float32
}

// EqualBase returns true if two input Float32Base structs, p and q, are 'equal', defined
// within a user-defined degree of precision.
func EqualBase(p Float32Base, q Float32Base, precision float32) bool {
	if !equalFloatPrecision(p.A, q.A, precision) {
		return false
	}
	if !equalFloatPrecision(p.C, q.C, precision) {
		return false
	}
	if !equalFloatPrecision(p.G, q.G, precision) {
		return false
	}
	if !equalFloatPrecision(p.T, q.T, precision) {
		return false
	}
	return true
}

// equalFloatPrecision is a helper function of EqualBase, which allows us to determine if
// two float32 are 'equal' within a user-specified degree of precision. Note that this is a
// relative precision (normalized to one of the values), and catches divide by zero errors.
func equalFloatPrecision(a float32, b float32, precision float32) bool {
	if a == 0 && b == 0 {
		return true
	}
	if a == 0 {
		return float32(math.Abs(float64(a-b)))/b < precision
	}
	return float32(math.Abs(float64(a-b)))/a < precision
}

// Entropy calculates the Shannon Entropy of a pDNA base.
func Entropy(base Float32Base) float64 {
	var answer float64 = 0.0

	// for each of these, we'll check p > 0 to avoid log2(0) calculation.
	if base.A > 0 {
		answer += -float64(base.A) * math.Log2(float64(base.A))
	}
	if base.C > 0 {
		answer += -float64(base.C) * math.Log2(float64(base.C))
	}
	if base.G > 0 {
		answer += -float64(base.G) * math.Log2(float64(base.G))
	}
	if base.T > 0 {
		answer += -float64(base.T) * math.Log2(float64(base.T))
	}

	return answer
}
