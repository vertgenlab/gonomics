package dnaThreeBit

import (
	"testing"
)

//func Append(fragment *ThreeBit, b ThreeBitBase) *ThreeBit { // too confusing to have Append and append (from std lib) in this package?
//func Cat(a *ThreeBit, b *ThreeBit) { // I was tempted to call this ligate
//func Copy(a *ThreeBit) *ThreeBit {

var dnaSeqs = []string{
	"TCATACGTTTTTTTTTTTTTCTGTC",
	"TCAAAACCCCCGGGGTTTTTCTGTC",
	"TCATACGTACGTACGTCCCCCTGCCCC",
	"TCATGGGGGGGGCCAGTACGTTGGCT",
	"TCATGGGGGGGGCCAGTACGTTGGCTTCAAAACCCCCGGGGTTTTTCTGTC",
}

func TestCopy(t *testing.T) {
	for _, currSeq := range dnaSeqs {
		original := FromString(currSeq)
		copyOf := Copy(original)
		original = Append(original, A) // mess with the original after it is copied
		stringAgain := ToString(copyOf)

		if stringAgain != currSeq {
			t.Errorf("Error: expected to get %s, but got %s instead.\n", currSeq, stringAgain)
		}
	}
}

var catTests = []struct {
	one  string // fragment one
	two  string // fragment two
	both string // what answer should be after cat
}{
	{"ACGT", "TTTGGCCAAA", "ACGTTTTGGCCAAA"},
	{"ACGT", "TTTGGCCAAACCCCTTTTGTGACGT", "ACGTTTTGGCCAAACCCCTTTTGTGACGT"},
	{"", "TTTGGCCAAA", "TTTGGCCAAA"},
	{"ACGT", "", "ACGT"},
}

func TestCat(t *testing.T) {
	for _, currTest := range catTests {
		tbOne := FromString(currTest.one)
		tbTwo := FromString(currTest.two)
		Cat(tbOne, tbTwo) // adds tbTwo to tbOne
		bothCalc := ToString(tbOne)

		if bothCalc != currTest.both {
			t.Errorf("Error: expected to get %s, but got %s instead.\n", currTest.both, bothCalc)
		}
	}
}
