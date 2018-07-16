package dna

import (
	"strings"
	"testing"
)

var equalStringsAndBases = []struct {
	characters string // first input
	bases      []Base // second input
}{
	{"ACGTNacgtn", []Base{A, C, G, T, N, a, c, g, t, n}},
	{"ACACGGGGNNNNNN", []Base{A, C, A, C, G, G, G, G, N, N, N, N, N, N}},
	{"NNNNnnnnTGTGTG", []Base{N, N, N, N, n, n, n, n, T, G, T, G, T, G}},
	{"NNNNgtgtacaNNN", []Base{N, N, N, N, g, t, g, t, a, c, a, N, N, N}},
	{"nnnacgttgcannnaacc", []Base{n, n, n, a, c, g, t, t, g, c, a, n, n, n, a, a, c, c}},
}

func TestStringToBases(t *testing.T) {
	for _, test := range equalStringsAndBases {
		actual, err := StringToBases(test.characters)
		if CompareSeqsCaseSensitive(actual, test.bases) != 0 || err != nil {
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
