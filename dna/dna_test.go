package dna

import "testing"

var dnaStrings = []string{
	"ACGTacgtNn",
	"AAAAAACCCCCGGGGTTTTT",
	"aaaaaCCCCTTTaaaaa",
	"NNNNNAAAaaaTTT",
}

func TestDnaToFromString(t *testing.T) {
	for _, input := range dnaStrings {
		bases := StringToBases(input)
		answer := BasesToString(bases)
		if input != answer {
			t.Errorf("Converting %s to bases and back gave %s", input, answer)
		}
	}
}

var rangeUpperTests = []struct {
	input    []Base // input
	start    int
	end      int
	expected []Base // after applying uppercase
}{
	{[]Base{A, C, G, T, N, a, c, g, t, n}, 6, 8, []Base{A, C, G, T, N, a, C, G, t, n}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, 6, 8, []Base{A, C, G, T, N, A, C, G, T, N}},
	{[]Base{a, c, g, t, n, A, C, G, T, N}, 0, 1, []Base{A, c, g, t, n, A, C, G, T, N}},
	{[]Base{a, c, g, t, n, a, c, g, t, n}, 7, 10, []Base{a, c, g, t, n, a, c, G, T, N}},
	{[]Base{a, a, a, g, g, t, t, t, n, N, C, g, t, n}, 1, 4, []Base{a, A, A, G, g, t, t, t, n, N, C, g, t, n}},
}

func TestRangeUpper(t *testing.T) {
	for _, test := range rangeUpperTests {
		actual := make([]Base, len(test.input))
		copy(actual, test.input)
		RangeToUpper(actual, test.start, test.end)
		if CompareSeqsCaseSensitive(actual, test.expected) != 0 {
			t.Errorf("Converting %v to uppercase from %d to %d gave %v when %v was expected.", test.input, test.start, test.end, actual, test.expected)
		}
	}
}

var rangeLowerTests = []struct {
	input    []Base // input
	start    int
	end      int
	expected []Base // after applying uppercase
}{
	{[]Base{A, C, G, T, N, a, c, g, t, n}, 6, 8, []Base{A, C, G, T, N, a, c, g, t, n}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, 6, 8, []Base{A, C, G, T, N, A, c, g, T, N}},
	{[]Base{A, C, g, t, n, A, C, G, T, N}, 0, 1, []Base{a, C, g, t, n, A, C, G, T, N}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, 7, 10, []Base{A, C, G, T, N, A, C, g, t, n}},
	{[]Base{a, a, A, G, G, t, t, t, n, N, C, g, t, n}, 1, 4, []Base{a, a, a, g, G, t, t, t, n, N, C, g, t, n}},
}

func TestRangeLower(t *testing.T) {
	for _, test := range rangeLowerTests {
		actual := make([]Base, len(test.input))
		copy(actual, test.input)
		RangeToLower(actual, test.start, test.end)
		if CompareSeqsCaseSensitive(actual, test.expected) != 0 {
			t.Errorf("Converting %v to lowercase from %d to %d gave %v when %v was expected.", test.input, test.start, test.end, actual, test.expected)
		}
	}
}

var allUpperTests = []struct {
	input    []Base // input
	expected []Base // after applying uppercase
}{
	{[]Base{A, C, G, T, N, a, c, g, t, n}, []Base{A, C, G, T, N, A, C, G, T, N}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, []Base{A, C, G, T, N, A, C, G, T, N}},
	{[]Base{a, c, g, t, n, A, C, G, T, N}, []Base{A, C, G, T, N, A, C, G, T, N}},
	{[]Base{a, c, g, t, n, a, c, g, t, n}, []Base{A, C, G, T, N, A, C, G, T, N}},
	{[]Base{a, a, a, g, g, t, t, t, n, N, C, g, t, n}, []Base{A, A, A, G, G, T, T, T, N, N, C, G, T, N}},
}

func TestAllUpper(t *testing.T) {
	for _, test := range allUpperTests {
		actual := make([]Base, len(test.input))
		copy(actual, test.input)
		AllToUpper(actual)
		if CompareSeqsCaseSensitive(actual, test.expected) != 0 {
			t.Errorf("Converting %v to uppercase gave %v when %v was expected.", test.input, actual, test.expected)
		}
	}
}

var allLowerTests = []struct {
	input    []Base // input
	expected []Base // after applying uppercase
}{
	{[]Base{A, C, G, T, N, a, c, g, t, n}, []Base{a, c, g, t, n, a, c, g, t, n}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, []Base{a, c, g, t, n, a, c, g, t, n}},
	{[]Base{a, c, g, t, n, A, C, G, T, N}, []Base{a, c, g, t, n, a, c, g, t, n}},
	{[]Base{a, c, g, t, n, a, c, g, t, n}, []Base{a, c, g, t, n, a, c, g, t, n}},
	{[]Base{A, a, a, g, g, T, T, T, n, n, c, g, t, N}, []Base{a, a, a, g, g, t, t, t, n, n, c, g, t, n}},
}

func TestAllLower(t *testing.T) {
	for _, test := range allLowerTests {
		actual := make([]Base, len(test.input))
		copy(actual, test.input)
		AllToLower(actual)
		if CompareSeqsCaseSensitive(actual, test.expected) != 0 {
			t.Errorf("Converting %v to lowercase gave %v when %v was expected.", test.input, actual, test.expected)
		}
	}
}

var comparisonTests = []struct {
	inputOne              []Base // first seq
	inputTwo              []Base // second seq
	expectedIgnoringCase  int    // expected result when ignoring case
	expectedCaseSensitive int    // expected result when considering case
}{
	{[]Base{A, C, G, T, N, a, c, g, t, n}, []Base{a, c, g, t, n, a, c, g, t, n}, 0, -1},
	{[]Base{a, C, G, T, N, a, c, g, t, n}, []Base{A, c, g, t, n, a, c, g, t, n}, 0, 1},
	{[]Base{T, c, g, t, n, A, C, G, T, N}, []Base{c, c, g, t, n, a, c, g, t, n}, 1, -1},
	{[]Base{a, c, g, t, n, a, c, g, t, n}, []Base{a, T, g, t, n, a, c, g, t, n}, -1, 1},
}

func TestComparison(t *testing.T) {
	for _, test := range comparisonTests {
		actual := CompareSeqsIgnoreCase(test.inputOne, test.inputTwo)
		if CompareSeqsIgnoreCase(test.inputOne, test.inputTwo) != test.expectedIgnoringCase {
			t.Errorf("The comparison of %v and %v when ignoring case gave %d when %d was expected", test.inputOne, test.inputTwo, actual, test.expectedIgnoringCase)
		}
		actual = CompareSeqsCaseSensitive(test.inputOne, test.inputTwo)
		if CompareSeqsCaseSensitive(test.inputOne, test.inputTwo) != test.expectedCaseSensitive {
			t.Errorf("The case-sensitive comparison of %v and %v gave %d when %d was expected", test.inputOne, test.inputTwo, actual, test.expectedCaseSensitive)
		}
	}
}
