package dnaThreeBit

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"testing"
)

//func NewThreeBitRainbow(inSeq []dna.Base, padding ThreeBitBase) []*ThreeBit {

var dnaFragments = []string{
	"TCATACGTTTTTTTTTTTTTCTGTC",
	"TCAAAACCCCCGGGGTTTTTCTGTC",
	"TCATACGTACGTACGTCCCCCTGCCCC",
	"TCATGGGGGGGGCCAGTACGTTGGCT",
	"TCATGGGGGGGGCCAGTACGTTGGCTTCAAAACCCCCGGGGTTTTTCTGTC",
}

func TestRainbow(t *testing.T) {
	for _, currString := range dnaFragments {
		currBases, err := dna.StringToBases(currString)
		if err != nil {
			log.Panicf("error converting string to Bases")
		}
		rainbowTable := NewThreeBitRainbow(currBases, G)

		for i, color := range rainbowTable {
			shouldBeOriginal := dna.BasesToString(RangeToDnaBases(color, i, i+len(currString)))
			if shouldBeOriginal != currString {
				t.Errorf("Error: expected to get %s, but got %s instead.\n", currString, shouldBeOriginal)
			}
		}
	}
}
