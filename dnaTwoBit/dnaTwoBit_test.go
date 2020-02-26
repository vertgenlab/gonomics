package dnaTwoBit

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

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
		frag := NewTwoBit(bases)
		singleBase = GetBase(frag, 0)
		if singleBase != dna.T {
			t.Errorf("Error: expected to get a T, but got a %c. %64b\n", dna.BaseToRune(singleBase), frag.Seq[0])
		}
		singleBase = GetBase(frag, 1)
		if singleBase != dna.C {
			t.Errorf("Error: expected to get a C, but got a %c. %64b\n", dna.BaseToRune(singleBase), frag.Seq[0])
		}
		singleBase = GetBase(frag, 2)
		if singleBase != dna.A {
			t.Errorf("Error: expected to get an A, but got a %c. %64b\n", dna.BaseToRune(singleBase), frag.Seq[0])
		}
		singleBase = GetBase(frag, 21)
		if singleBase != dna.T {
			t.Errorf("Error: expected to get a T, but got a %c. %64b\n", dna.BaseToRune(singleBase), frag.Seq[0])
		}
		singleBase = GetBase(frag, 24)
		if singleBase != dna.C {
			t.Errorf("Error: expected to get a C, but got a %c. %64b\n", dna.BaseToRune(singleBase), frag.Seq[0])
		}
	}
}
