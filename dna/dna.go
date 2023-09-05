// Package dna implements a data structure for storage and manipulation of sequences of DNA.
package dna

// Base stores a single nucleotide as a byte.
type Base byte

const (
	A      Base = 0
	C      Base = 1
	G      Base = 2
	T      Base = 3
	N      Base = 4
	LowerA Base = 5
	LowerC Base = 6
	LowerG Base = 7
	LowerT Base = 8
	LowerN Base = 9
	Gap    Base = 10
	Dot    Base = 11
	Nil    Base = 12
)

// CreatAllGaps creates a DNA sequence of Gap with length of numGaps.
func CreateAllGaps(numGaps int) []Base {
	answer := make([]Base, numGaps)
	for i := range answer {
		answer[i] = Gap
	}
	return answer
}

// CreatAllN creates a DNA sequence of N with length of numGaps.
func CreateAllNs(numGaps int) []Base {
	answer := make([]Base, numGaps)
	for i := range answer {
		answer[i] = N
	}
	return answer
}
