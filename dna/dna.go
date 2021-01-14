// Package dna provides functionality to create and manipulate DNA sequences.
package dna

// Base stores a single nucleotide as a byte.
type Base byte

const (
	A      Base = 'A'
	C      Base = 'C'
	G      Base = 'G'
	T      Base = 'T'
	N      Base = 'N'
	LowerA Base = 'a'
	LowerC Base = 'c'
	LowerG Base = 'g'
	LowerT Base = 't'
	LowerN Base = 'n'
	Gap    Base = '-'
	Dot    Base = '.'
	Nil    Base = '*'
)

// CreatAllGaps creates a DNA sequence of Gap with length of numGaps
func CreateAllGaps(numGaps int64) []Base {
	answer := make([]Base, numGaps)
	for i := range answer {
		answer[i] = Gap
	}
	return answer
}

// CreatAllN creates a DNA sequence of N with length of numGaps
func CreateAllNs(numGaps int64) []Base {
	answer := make([]Base, numGaps)
	for i := range answer {
		answer[i] = N
	}
	return answer
}
