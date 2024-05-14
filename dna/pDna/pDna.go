package pDna

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"math"
)

// Float32Base encodes a DNA base as a probability vector using float32 precision.
// Note that gap probabilities can be stored implicitly as 1 - (A + C + G + T).
type Float32Base struct {
	A float32
	C float32
	G float32
	T float32
}

/*
// Float32Diff encodes a probability vector using float64 precision that describes the difference between 2 pDNA probability vectors
// The probability of A, C, G, and T should add up to 0
type Float64Diff struct {
	A float64
	C float64
	G float64
	T float64
}
*/

// IsGap returns true if the DNA base's probability vector is 0 at all 4 bases, indicating that the base is a gap
func IsGap(p Float32Base) bool {
	if p.A == 0 && p.C == 0 && p.G == 0 && p.T == 0 {
		return true
	} else {
		return false
	}
}

/*
// IsNonGap returns true if the DNA base's probability vector adds up to 1, indicating that the base is a valid non-gap
func IsNonGap(p Float32Base) bool {
	if p.A+p.C+p.G+p.T == 1 {
		return true
	} else {
		return false
	}
}

// IsN returns true if the DNA base's probability vector is the same but not 0 (==0.25) at all 4 bases, indicating that the base is an N
func IsN(p Float32Base) bool {
	if p.A != 0 && p.A == p.C && p.A == p.G && p.A == p.T {
		return true
	} else {
		return false
	}
}
*/

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

/*
// Diff returns the vector that is the difference between 2 pDNA bases' probability vectors
func Diff(p Float32Base, q Float32Base) Float64Diff {
	return Float64Diff{
		A: float64(p.A - q.A),
		C: float64(p.C - q.C),
		G: float64(p.G - q.G),
		T: float64(p.T - q.T),
	}
}

// Mag returns the magnitude of a difference probability vector
// Mag is a float64
func Mag(d Float64Diff) float64 {
	return math.Sqrt(math.Pow(d.A, 2) + math.Pow(d.C, 2) + math.Pow(d.G, 2) + math.Pow(d.T, 2))
}

// Dist returns a score for the distance that separates 2 pDNA bases
// The distance score is the magnitude of the vector that is the difference between the 2 pDNA bases' probability vectors
// The distance score is a float64
func Dist(p Float32Base, q Float32Base) float64 {
	return Mag(Diff(p, q))
}
*/

// Dot
func Dot(p Float32Base, q Float32Base) float64 {
	return float64(p.A*q.A + p.C*q.C + p.G*q.G + p.T*q.T)
}

// DotSubstProb
// Make sure input is probabilities not likelihoods
func DotSubstProb(p Float32Base, q Float32Base) float64 {
	return 1 - Dot(p, q)
}

/*
// MakeValid converts the likelihoods of each base at a position of a sequence into a pdna for that position
func MakeValid(p Float32Base) Float32Base {
	total := p.A + p.C + p.G + p.T
	if total == 0 {
		return Float32Base{
			A: 0,
			C: 0,
			G: 0,
			T: 0,
		}
	} else {
		return Float32Base{
			A: p.A / total,
			C: p.C / total,
			G: p.G / total,
			T: p.T / total,
		}
	}
}
*/

func DnaToPdna(base dna.Base) Float32Base {
	var pbase Float32Base
	switch base {
	case 0: // A
		pbase.A = 1
		pbase.C = 0
		pbase.G = 0
		pbase.T = 0
	case 1: // C
		pbase.A = 0
		pbase.C = 1
		pbase.G = 0
		pbase.T = 0
	case 2: // G
		pbase.A = 0
		pbase.C = 0
		pbase.G = 1
		pbase.T = 0
	case 3: // T
		pbase.A = 0
		pbase.C = 0
		pbase.G = 0
		pbase.T = 1
	case 4: // N
		pbase.A = 0.25
		pbase.C = 0.25
		pbase.G = 0.25
		pbase.T = 0.25
	case 10: // Gap
		pbase.A = 0
		pbase.C = 0
		pbase.G = 0
		pbase.T = 0
	default:
		log.Fatalf("pDna.DnaToPdna encountered an inconvertible base.\n")
	}
	return pbase
}
