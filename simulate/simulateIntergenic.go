package simulate

import (
	"github.com/vertgenlab/gonomics/dna"
)

//RandIntergenicSeq makes a randomly generated DNA sequence of a specified length and GC content. Unlike RandGene, it does not have to be divisible by 3.
func RandIntergenicSeq(GcContent float64, lenSeq int) []dna.Base {
	var answer []dna.Base = make([]dna.Base, lenSeq)
	for i := range answer {
		answer[i] = chooseRandomBase(GcContent)
	}
	return answer
}

