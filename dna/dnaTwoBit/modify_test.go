package dnaTwoBit

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

var dnaSeqs = [][]dna.Base{
	dna.StringToBases("TCATACGTTTTTTTTTTTTTCTGTC"),
	dna.StringToBases("TCAAAACCCCCGGGGTTTTTCTGTC"),
	dna.StringToBases("TCATACGTACGTACGTCCCCCTGCCCC"),
	dna.StringToBases("TCATGGGGGGGGCCAGTACGTTGGCT"),
	dna.StringToBases("TCATGGGGGGGGCCAGTACGTTGGCTTCAAAACCCCCGGGGTTTTTCTGTC"),
}

func TestCopy(t *testing.T) {
	for _, currSeq := range dnaSeqs {
		original := NewTwoBit(currSeq)
		copyOf := Copy(original)
		original = Append(original, 0)
		stringAgain := ToString(copyOf)
		result := dna.BasesToString(currSeq)
		if ToString(original) != dna.BasesToString(append(currSeq, 0)) {
			t.Errorf("Error: expected to get %s, but got %s instead.\n", ToString(original), dna.BasesToString(append(currSeq, 0)))
		}
		if stringAgain != result {
			t.Errorf("Error: expected to get %s, but got %s instead.\n", result, stringAgain)
		}
	}
}

var catTests = []struct {
	one  []dna.Base // fragment one
	two  []dna.Base // fragment two
	both []dna.Base // what answer should be after cat
}{
	{dna.StringToBases("ACGT"), dna.StringToBases("TTTGGCCAAA"), dna.StringToBases("ACGTTTTGGCCAAA")},
	{dna.StringToBases("ACGT"), dna.StringToBases("TTTGGCCAAACCCCTTTTGTGACGT"), dna.StringToBases("ACGTTTTGGCCAAACCCCTTTTGTGACGT")},
	{[]dna.Base{}, dna.StringToBases("TTTGGCCAAA"), dna.StringToBases("TTTGGCCAAA")},
	{dna.StringToBases("ACGT"), []dna.Base{}, dna.StringToBases("ACGT")},
}

func TestCat(t *testing.T) {
	for _, currTest := range catTests {
		tbOne := NewTwoBit(currTest.one)
		tbTwo := NewTwoBit(currTest.two)
		Cat(tbOne, tbTwo) // adds tbTwo to tbOne
		bothCalc := ToString(tbOne)
		result := dna.BasesToString(currTest.both)
		if bothCalc != result {
			t.Errorf("Error: expected to get %s, but got %s instead.\n", currTest.both, result)
		}
	}
}
