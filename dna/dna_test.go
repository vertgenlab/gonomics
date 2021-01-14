package dna

import (
	"strings"
	"testing"
)

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
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, 6, 8, []Base{A, C, G, T, N, LowerA, C, G, LowerT, LowerN}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, 6, 8, []Base{A, C, G, T, N, A, C, G, T, N}},
	{[]Base{LowerA, LowerC, LowerG, LowerT, LowerN, A, C, G, T, N}, 0, 1, []Base{A, LowerC, LowerG, LowerT, LowerN, A, C, G, T, N}},
	{[]Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}, 7, 10, []Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, G, T, N}},
	{[]Base{LowerA, LowerA, LowerA, LowerG, LowerG, LowerT, LowerT, LowerT, LowerN, N, C, LowerG, LowerT, LowerN}, 1, 4, []Base{LowerA, A, A, G, LowerG, LowerT, LowerT, LowerT, LowerN, N, C, LowerG, LowerT, LowerN}},
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
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, 6, 8, []Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, 6, 8, []Base{A, C, G, T, N, A, LowerC, LowerG, T, N}},
	{[]Base{A, C, LowerG, LowerT, LowerN, A, C, G, T, N}, 0, 1, []Base{LowerA, C, LowerG, LowerT, LowerN, A, C, G, T, N}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, 7, 10, []Base{A, C, G, T, N, A, C, LowerG, LowerT, LowerN}},
	{[]Base{LowerA, LowerA, A, G, G, LowerT, LowerT, LowerT, LowerN, N, C, LowerG, LowerT, LowerN}, 1, 4, []Base{LowerA, LowerA, LowerA, LowerG, G, LowerT, LowerT, LowerT, LowerN, N, C, LowerG, LowerT, LowerN}},
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
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, []Base{A, C, G, T, N, A, C, G, T, N}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, []Base{A, C, G, T, N, A, C, G, T, N}},
	{[]Base{LowerA, LowerC, LowerG, LowerT, LowerN, A, C, G, T, N}, []Base{A, C, G, T, N, A, C, G, T, N}},
	{[]Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}, []Base{A, C, G, T, N, A, C, G, T, N}},
	{[]Base{LowerA, LowerA, LowerA, LowerG, LowerG, LowerT, LowerT, LowerT, LowerN, N, C, LowerG, LowerT, LowerN}, []Base{A, A, A, G, G, T, T, T, N, N, C, G, T, N}},
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
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, []Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, []Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}},
	{[]Base{LowerA, LowerC, LowerG, LowerT, LowerN, A, C, G, T, N}, []Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}},
	{[]Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}, []Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}},
	{[]Base{A, LowerA, LowerA, LowerG, LowerG, T, T, T, LowerN, LowerN, LowerC, LowerG, LowerT, N}, []Base{LowerA, LowerA, LowerA, LowerG, LowerG, LowerT, LowerT, LowerT, LowerN, LowerN, LowerC, LowerG, LowerT, LowerN}},
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
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, []Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}, 0, -1},
	{[]Base{LowerA, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, []Base{A, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}, 0, 1},
	{[]Base{T, LowerC, LowerG, LowerT, LowerN, A, C, G, T, N}, []Base{LowerC, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}, 1, -1},
	{[]Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}, []Base{LowerA, T, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}, -1, 1},
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
