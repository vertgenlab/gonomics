package pDna

import (
	"math"
	"math/rand"
)

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

// Scale multiplies the four values in a pDNA base (A, C, T, G) by the multiplier
func Scale(pdnaBase Float32Base, multiplier float32) Float32Base {
	return Float32Base{A: pdnaBase.A*multiplier, C: pdnaBase.C*multiplier, G: pdnaBase.G*multiplier, T: pdnaBase.T*multiplier}
}

// Sum adds the respective four values in two pDNA bases (A, C, T, G). Permits bases with values greater than 1. 
func Sum(base1 Float32Base, base2 Float32Base) Float32Base {
	return Float32Base{A: base1.A+base2.A, C: base1.C+base2.C, G: base1.G+base2.G, T: base1.T+base2.T}
}

// SumsToOne checks if the total probabilities of a pDNA base sum to 1 
func SumsToOne(base Float32Base, precision float32) bool {
	if !equalFloatPrecision(base.A + base.C + base.G + base.T, 1, precision) {
		return false
	}
	return true
}

// Random randomly generates a pDNA base that sums to 1
func RandBase(seedSet bool, setSeed int64) Float32Base {
	var answer Float32Base
	if !seedSet {
		rand.Seed(setSeed)
	} 
	sumsToOne := false
	for !sumsToOne {
		aFloat := rand.Float32()
		cFloat := rand.Float32()
		gFloat := rand.Float32()
		tFloat := rand.Float32()
		sum := aFloat + cFloat + gFloat + tFloat
		answer.A = aFloat / sum
		answer.C = cFloat / sum
		answer.G = gFloat / sum
		answer.T = tFloat / sum
		if SumsToOne(answer, 1e-3) {
			sumsToOne = true
			break
		} 
	}
	
	
	return answer
}
