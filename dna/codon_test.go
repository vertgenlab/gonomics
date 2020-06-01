package dna

import "testing"

var expected []*Codon = []*Codon{{StringToBases("AAG")}, {StringToBases("TGA")}, {StringToBases("AAC")}}
var testSeq []Base = StringToBases("AAGTGAAAC")

var test2DNA []Base = StringToBases("atgcagatttttgtgaaaaccctgaccggcaaaacctga")
var expectedShortProtein string = "MQIFVKTLTGKT*"
var expectedProtein string = "MetGlnIlePheValLysThrLeuThrGlyLysThrTer"

func TestBasesToCodon(t *testing.T) {
	input := BasesToCodons(testSeq)
	for i := 0; i < len(input); i++ {
		if BasesToString(input[i].Seq) != BasesToString(expected[i].Seq) {
			t.Errorf("Do not match. Input: %s. Expected: %s.", BasesToString(input[i].Seq), BasesToString(expected[i].Seq))
		}
	}
}

func TestTranslateToString(t *testing.T) {
	input := TranslateToString(test2DNA)
	if input != expectedProtein {
		t.Errorf("Do not match. Input : %s. Expected: %s.", input, expectedProtein)
	}
}

func TestTranslateToShortString(t *testing.T) {
	input := TranslateToShortString(test2DNA)
	if input != expectedShortProtein {
		t.Errorf("Do not match. Input : %s. Expected: %s.", input, expectedShortProtein)
	}
}
