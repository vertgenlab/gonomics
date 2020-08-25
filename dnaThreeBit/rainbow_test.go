package dnaThreeBit

import (
	"github.com/vertgenlab/gonomics/dna"
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
		rainbowTable := NewThreeBitRainbow(dna.StringToBases(currString), G)

		for i, color := range rainbowTable {
			shouldBeOriginal := dna.BasesToString(SectionToDnaBases(color, i, i+len(currString)))
			if shouldBeOriginal != currString {
				t.Errorf("Error: expected to get %s, but got %s instead.\n", currString, shouldBeOriginal)
			}
		}
	}
}
