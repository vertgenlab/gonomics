package fasta

import (
	"testing"
	"github.com/vertgenlab/gonomics/dna"
	"sort"
)

var contigOneSeq = dna.StringToBases("ACGTGAGTGAGTAGGACCACGATGACACGANNTGA")
var contigTwoSeq = dna.StringToBases("GgtAC")
var contigThreeSeq = dna.StringToBases("GTAGTGAGTGA")
var inputFasta []*Fasta = []*Fasta{{"apple", contigOneSeq}, {"banana", contigTwoSeq}, {"carrot", contigThreeSeq}}
var expected []int = []int{1, 2, 3, 11, 30}

func TestMakeContigList(t *testing.T) {
	input := makeContigList(inputFasta, true)
	sort.Ints(input)
	if !ContigAllAreEqual(input, expected) {
			t.Errorf("Do not match. Input: %v. Expected: %v. \n", input, expected)
		}
}

func ContigAllAreEqual(input []int, expected []int) bool {
	for i := 0; i < len(input); i++ {
		if input[i] != expected[i] {
			return false
		}
	}
	return true
}