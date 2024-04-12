package pDna

import (
	"math"
)

// Float32Base encodes a DNA base as a probability vector using float32 precision.
// Note that gap probabilities can be stored implicitly as 1 - (A + C + G + T).
// The probabilities of A, C, G, and T should add up to 1
type Float32Base struct {
	A float32
	C float32
	G float32
	T float32
}

// Float32Diff encodes a probability vector using float64 precision that describes the difference between 2 pDNA probability vectors
// The probability of A, C, G, and T should add up to 0
type Float64Diff struct {
	A float64
	C float64
	G float64
	T float64
}

// IsGap returns true if the DNA base's probability vector is 0 at all 4 bases, indicating that the base is a gap
func IsGap(p Float32Base) bool {
	if p.A == 0 && p.C == 0 && p.G == 0 && p.T == 0 {
		//fmt.Printf("found gap. p: %v\n", p) // TODO: remove after debugging. Check isGap works
		return true
	} else {
		return false
	}
}

// TODO: log fatal instead of return false here? Update pFasta.go QC accordingly
// IsNonGap returns true if the DNA base's probability vector adds up to 1, indicating that the base is a valid non-gap
func IsNonGap(p Float32Base) bool {
	if p.A+p.C+p.G+p.T == 1 {
		return true
	} else {
		return false
	}
}

// TODO: should I worry about precision here? Use EqualBase, with some standard precision?
// IsN returns true if the DNA base's probability vector is the same but not 0 (==0.25) at all 4 bases, indicating that the base is an N
func IsN(p Float32Base) bool {
	if p.A != 0 && p.A == p.C && p.A == p.G && p.A == p.T {
		//fmt.Printf("found N. p: %v\n", p) // TODO: remove after debugging. Check isN works
		return true
	} else {
		return false
	}
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
	//return Mag(Diff(p, q)) // regular version

	// logspace version
	// p and q are log10() numbers
	// convert to regular numbers with math.Pow
	pUnlog := []float64{
		math.Pow(10, float64(p.A)),
		math.Pow(10, float64(p.C)),
		math.Pow(10, float64(p.G)),
		math.Pow(10, float64(p.T)),
	}
	qUnlog := []float64{
		math.Pow(10, float64(q.A)),
		math.Pow(10, float64(q.C)),
		math.Pow(10, float64(q.G)),
		math.Pow(10, float64(q.T)),
	}
	// calculate diff
	diff := []float64{
		pUnlog[0] - qUnlog[0],
		pUnlog[1] - qUnlog[1],
		pUnlog[2] - qUnlog[2],
		pUnlog[3] - qUnlog[3],
	}
	// calculate mag
	mag := math.Sqrt(math.Pow(diff[0], 2) + math.Pow(diff[1], 2) + math.Pow(diff[2], 2) + math.Pow(diff[3], 2))
	return mag
}

// Dot
func Dot(p Float32Base, q Float32Base) float64 {
	// TODO: remove the below if loop after debugging? Should have already checked in pfaFindFast/efficient.go for gap, so this is not gap
	//if (p.A+p.C+p.G+p.T != 1) || (q.A+q.C+q.G+q.T != 1) {
	//	log.Fatalf("p or q sum not 1. p:%v, sum: %v, q: %v, sum: %v\n", p, p.A+p.C+p.G+p.T, q, q.A+q.C+q.G+q.T)
	//}
	return float64(p.A*q.A + p.C*q.C + p.G*q.G + p.T*q.T)
}

// DotSubstProb
// Make sure input is probabilities not likelihoods
func DotSubstProb(p Float32Base, q Float32Base) float64 {
	//return 1 - Dot(p, q) // regular version

	// logspace version
	// p and q are log10() numbers
	// convert to regular numbers with math.Pow
	pUnlog := []float64{
		math.Pow(10, float64(p.A)),
		math.Pow(10, float64(p.C)),
		math.Pow(10, float64(p.G)),
		math.Pow(10, float64(p.T)),
	}
	qUnlog := []float64{
		math.Pow(10, float64(q.A)),
		math.Pow(10, float64(q.C)),
		math.Pow(10, float64(q.G)),
		math.Pow(10, float64(q.T)),
	}
	dot := pUnlog[0]*qUnlog[0] + pUnlog[1]*qUnlog[1] + pUnlog[2]*qUnlog[2] + pUnlog[3]*qUnlog[3]
	return 1 - dot
}

// TODO: maybe replace LikelihoodsToPdna in reconstruct/reconstruct.go in branch reconstructSeq_nodeMods
// MakeValid converts the likelihoods of each base at a position of a sequence into a pdna for that position
func MakeValid(p Float32Base) Float32Base {
	total := p.A + p.C + p.G + p.T
	//fmt.Printf("raw ACGT: %g, %g, %g, %g\n", p.A, p.C, p.G, p.T)
	//fmt.Printf("raw total: %g\n", total)
	//fmt.Printf("normalized ACGT: %g, %g, %g, %g\n", p.A/total, p.C/total, p.G/total, p.T/total)
	//fmt.Printf("normalized total: %g\n", p.A/total+p.C/total+p.G/total+p.T/total)
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
