package dnaThreeBit

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

//func RuneToThreeBitBase(r rune) ThreeBitBase {
//func ThreeBitBaseToRune(base ThreeBitBase) rune {
//func ThreeBitBaseToString(b ThreeBitBase) string {
//func FromString(s string) *ThreeBit {
//func ToString(fragment *ThreeBit) string {
//func SectionToDnaBases(fragment *ThreeBit, start int, end int) []dna.Base {
//func ToDnaBases(fragment *ThreeBit) []dna.Base {

var dnaSamples = []string{
	"TCATACGTTTTTTTTTTTTTCTGTC",
	"TCAAAACCCCCGGGGTTTTTCTGTC",
	"TCATACGTACGTACGTCCCCCTGCCCC",
	"TCATGGGGGGGGCCAGTACGTTGGCT",
	"TCATGGGGGGGGCCAGTACGTTGGCTTCAAAACCCCCGGGGTTTTTCTGTC",
}

func TestFromStringAndBack(t *testing.T) {
	for _, input := range dnaSamples {
		threeBitVersion := FromString(input)
		stringAgain := ToString(threeBitVersion)

		if stringAgain != input {
			t.Errorf("Error: expected to get %s, but got %s instead.\n", input, stringAgain)
		}
	}
}

func TestFromStringAndBackViaDnaBase(t *testing.T) {
	for _, input := range dnaStrings {
		threeBitVersion := FromString(input)
		stringAgain := dna.BasesToString(ToDnaBases(threeBitVersion))

		if stringAgain != input {
			t.Errorf("Error: expected to get %s, but got %s instead.\n", input, stringAgain)
		}
	}
}

func TestIndividualBasesViaThreeBitBase(t *testing.T) {
	for _, input := range dnaStrings {
		tripleBit := FromString(input)
		singleBase := ThreeBitBaseToString(GetThreeBitBase(tripleBit, 0))
		if singleBase != "T" {
			t.Errorf("Error: expected to get a T, but got a %s. %64b\n", singleBase, tripleBit.Seq[0])
		}
		singleBase = ThreeBitBaseToString(GetThreeBitBase(tripleBit, 1))
		if singleBase != "C" {
			t.Errorf("Error: expected to get a C, but got a %s. %64b\n", singleBase, tripleBit.Seq[0])
		}
		singleBase = ThreeBitBaseToString(GetThreeBitBase(tripleBit, 2))
		if singleBase != "A" {
			t.Errorf("Error: expected to get an A, but got a %s. %64b\n", singleBase, tripleBit.Seq[0])
		}
		singleBase = ThreeBitBaseToString(GetThreeBitBase(tripleBit, 21))
		if singleBase != "T" {
			t.Errorf("Error: expected to get a T, but got a %s. %64b\n", singleBase, tripleBit.Seq[0])
		}
		singleBase = ThreeBitBaseToString(GetThreeBitBase(tripleBit, 24))
		if singleBase != "C" {
			t.Errorf("Error: expected to get a C, but got a %s. %64b\n", singleBase, tripleBit.Seq[0])
		}
	}
}
