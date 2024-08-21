package fasta

import (
	"fmt"
	"sort"
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

var contigOneSeq = dna.StringToBases("ACGTGAGTGAGTAGGACCACGATGACACGANNTGA")
var contigTwoSeq = dna.StringToBases("GgtAC")
var contigThreeSeq = dna.StringToBases("GTAGTGAGTGA")
var inputFasta []Fasta = []Fasta{{"apple", contigOneSeq}, {"banana", contigTwoSeq}, {"carrot", contigThreeSeq}}
var expected []int = []int{1, 2, 3, 11, 30}

func TestMakeContigList(t *testing.T) {
	input := MakeContigList(inputFasta, true)
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

func TestCalculateGenomeLength(t *testing.T) {
	input := MakeContigList(inputFasta, false)
	if calculateGenomeLength(input) != 49 {
		t.Errorf("Error in calculateGenomeLength. Expected: 49. Observed: %d", calculateGenomeLength(input))
	}
}

func TestCalculateN50L50(t *testing.T) {
	exp := []int{30, 1}
	input := MakeContigList(inputFasta, true)
	sort.Ints(input)
	fmt.Println(input)
	N50, L50 := CalculateN50L50(input, calculateGenomeLength(input)/2)
	if N50 != exp[0] || L50 != exp[1] {
		t.Errorf("Error in CalculateN50L50: Expected N50: %d L50: %d. Observed N50: %d L50: %d", exp[0], exp[1], N50, L50)
	}

}
