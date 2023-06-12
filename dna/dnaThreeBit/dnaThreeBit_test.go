package dnaThreeBit

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

//func BasesToUint64(seq []dna.Base, start int, end int, padding ThreeBitBase) uint64 {
//func basesToUint64WithOffset(seq []dna.Base, start int, end int, padding ThreeBitBase, offset int) uint64 {
//func GetThreeBitBase(fragment *ThreeBit, pos int) ThreeBitBase {
//func GetBase(fragment *ThreeBit, pos int) dna.Base {
//func NewThreeBit(inSeq []dna.Base, padding ThreeBitBase) *ThreeBit {
//func newThreeBitWithOffset(inSeq []dna.Base, padding ThreeBitBase, offset int) *ThreeBit {

var dnaStrings = []string{
	"TCATACGTTTTTTTTTTTTTCTGTC",
	"TCAAAACCCCCGGGGTTTTTCTGTC",
	"TCATACGTACGTACGTCCCCCTGCCCC",
	"TCATGGGGGGGGCCAGTACGTTGGCT",
}

func TestGetBaseWithOffset(t *testing.T) {
	var singleBase dna.Base
	for _, input := range dnaStrings {
		bases := dna.StringToBases(input)
		offset := 1
		tripleBit := newThreeBitWithOffset(bases, PaddingOne, offset)
		singleBase = GetBase(tripleBit, 0+offset)
		if singleBase != dna.T {
			t.Errorf("Error: expected to get a T, but got a %c. %64b\n", dna.BaseToRune(singleBase), tripleBit.Seq[0])
		}
		singleBase = GetBase(tripleBit, 1+offset)
		if singleBase != dna.C {
			t.Errorf("Error: expected to get a C, but got a %c. %64b\n", dna.BaseToRune(singleBase), tripleBit.Seq[0])
		}
		singleBase = GetBase(tripleBit, 2+offset)
		if singleBase != dna.A {
			t.Errorf("Error: expected to get an A, but got a %c. %64b\n", dna.BaseToRune(singleBase), tripleBit.Seq[0])
		}
		singleBase = GetBase(tripleBit, 21+offset)
		if singleBase != dna.T {
			t.Errorf("Error: expected to get a T, but got a %c. %64b\n", dna.BaseToRune(singleBase), tripleBit.Seq[0])
		}
		singleBase = GetBase(tripleBit, 24+offset)
		if singleBase != dna.C {
			t.Errorf("Error: expected to get a C, but got a %c. %64b\n", dna.BaseToRune(singleBase), tripleBit.Seq[0])
		}
	}
}

func TestGetBase(t *testing.T) {
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
