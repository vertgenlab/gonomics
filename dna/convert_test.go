package dna

import (
	"strings"
	"testing"
)

var equalStringsAndBases = []struct {
	characters string // first input
	bases      []Base // second input
}{
	{"ACGTNacgtn", []Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}},
	{"ACACGGGGNNNNNN", []Base{A, C, A, C, G, G, G, G, N, N, N, N, N, N}},
	{"NNNNnnnnTGTGTG", []Base{N, N, N, N, LowerN, LowerN, LowerN, LowerN, T, G, T, G, T, G}},
	{"NNNNgtgtacaNNN", []Base{N, N, N, N, LowerG, LowerT, LowerG, LowerT, LowerA, LowerC, LowerA, N, N, N}},
	{"nnnacgttgcannnaacc", []Base{LowerN, LowerN, LowerN, LowerA, LowerC, LowerG, LowerT, LowerT, LowerG, LowerC, LowerA, LowerN, LowerN, LowerN, LowerA, LowerA, LowerC, LowerC}},
}

func TestStringToBases(t *testing.T) {
	for _, test := range equalStringsAndBases {
		actual := StringToBases(test.characters)
		if CompareSeqsCaseSensitive(actual, test.bases) != 0 {
			t.Errorf("StringToBases(%s): expected %v, actual %v", test.characters, test.bases, actual)
		}
	}
}

func TestBasesToString(t *testing.T) {
	for _, test := range equalStringsAndBases {
		actual := BasesToString(test.bases)
		if strings.Compare(actual, test.characters) != 0 {
			t.Errorf("BasesToString(%v): expected %s, actual %v", test.bases, test.characters, actual)
		}
	}
}
