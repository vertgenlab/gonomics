package dna

import ()

type Base byte

// TODO: change these names so that all variables can be seen
// from outside.  Maybe BaseA, Basea, etc
const (
	A   Base = 0
	C   Base = 1
	G   Base = 2
	T   Base = 3
	N   Base = 4
	a   Base = 5
	c   Base = 6
	g   Base = 7
	t   Base = 8
	n   Base = 9
	Gap Base = 10
	Dot Base = 11
)

func CreateAllGaps(numGaps int64) []Base {
	answer := make([]Base, numGaps)
	for i := 0; i < len(answer); i++ {
		answer[i] = Gap
	}
	return answer
}
