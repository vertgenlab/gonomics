package dna

import "testing"

var revCompTests = []struct {
	input    []Base // input
	expected []Base // after applying reverseComplement
}{
	{[]Base{A, C, G, T, N, a, c, g, t, n}, []Base{n, a, c, g, t, N, A, C, G, T}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, []Base{N, A, C, G, T, N, A, C, G, T}},
	{[]Base{a, c, g, t, n, A, C, G, T, N}, []Base{N, A, C, G, T, n, a, c, g, t}},
	{[]Base{a, c, g, t, n, a, c, g, t, n}, []Base{n, a, c, g, t, n, a, c, g, t}},
	{[]Base{A, A, A, G, G, T, T, T, N, n, c, g, t, n}, []Base{n, a, c, g, n, N, A, A, A, C, C, T, T, T}},
}

func TestRevComp(t *testing.T) {
	for _, test := range revCompTests {
		actual := make([]Base, len(test.input))
		copy(actual, test.input)
		ReverseComplement(actual)
		if CompareSeqsCaseSensitive(actual, test.expected) != 0 {
			t.Errorf("Reverse complementing %v gave %v when %v was expected.", test.input, actual, test.expected)
		}
	}
}
