package dnaThreeBit

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

//func GetThreeBitBase(frag *ThreeBit, pos uint) dna.Base {
//func NewThreeBit(inSeq []dna.Base) *ThreeBit {

var dnaStrings = []string{
	"TCATACGTTTTTTTTTTTTTCTGTC",
	"TCAAAACCCCCGGGGTTTTTCTGTC",
	"TCATACGTACGTACGTCCCCCTGCCCC",
	"TCATGGGGGGGGCCAGTACGTTGGCT",
}

func TestDnaToFromString(t *testing.T) {
	var singleBase dna.Base
	for _, input := range dnaStrings {
		bases := dna.StringToBases(input)
		tripleBit := NewThreeBit(bases, PaddingOne)
		singleBase = GetBase(tripleBit, 0)
		if singleBase != dna.T {
			t.Errorf("Error: expected to get a T, but got a %c. %64b\n", dna.BaseToRune(singleBase), tripleBit.Seq[0])
		}
		singleBase = GetBase(tripleBit, 1)
		if singleBase != dna.C {
			t.Errorf("Error: expected to get a C, but got a %c. %64b\n", dna.BaseToRune(singleBase), tripleBit.Seq[0])
		}
		singleBase = GetBase(tripleBit, 2)
		if singleBase != dna.A {
			t.Errorf("Error: expected to get an A, but got a %c. %64b\n", dna.BaseToRune(singleBase), tripleBit.Seq[0])
		}
		singleBase = GetBase(tripleBit, 21)
		if singleBase != dna.T {
			t.Errorf("Error: expected to get a T, but got a %c. %64b\n", dna.BaseToRune(singleBase), tripleBit.Seq[0])
		}
		singleBase = GetBase(tripleBit, 24)
		if singleBase != dna.C {
			t.Errorf("Error: expected to get a C, but got a %c. %64b\n", dna.BaseToRune(singleBase), tripleBit.Seq[0])
		}
	}
}
