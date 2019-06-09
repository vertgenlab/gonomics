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

var deletionTests = []struct {
	input    []Base // input
	delStart int64
	delEnd   int64
	expected []Base // after applying deletion
}{
	{[]Base{A, C, G, T, N, a, c, g, t, n}, 0, 2, []Base{G, T, N, a, c, g, t, n}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, 7, 10, []Base{A, C, G, T, N, A, C}},
	{[]Base{a, c, g, t, n, A, C, G, T, N}, 3, 4, []Base{a, c, g, n, A, C, G, T, N}},
	{[]Base{a, c, g, t, n, a, c, g, t, n}, 4, 8, []Base{a, c, g, t, t, n}},
}

var insertionTests = []struct {
	seq      []Base // input
	pos      int64
	insSeq   []Base
	expected []Base // after applying deletion
}{
	{[]Base{A, C, G, T, N, a, c, g, t, n}, 0, []Base{T, A, T, A}, []Base{T, A, T, A, A, C, G, T, N, a, c, g, t, n}},
	{[]Base{A, C, G, T, N, a, c, g, t, n}, 10, []Base{T, A, T, A}, []Base{A, C, G, T, N, a, c, g, t, n, T, A, T, A}},
	{[]Base{A, C, G, T, N, a, c, g, t, n}, 3, []Base{A, C, G, T}, []Base{A, C, G, A, C, G, T, T, N, a, c, g, t, n}},
	{[]Base{A, C, G, T, N, a, c, g, t, n}, 4, []Base{A, C, G, T}, []Base{A, C, G, T, A, C, G, T, N, a, c, g, t, n}},
}

var replaceTests = []struct {
	seq      []Base // input
	start    int64
	end      int64
	insSeq   []Base
	expected []Base // after applying deletion
}{
	{[]Base{A, C, G, T, N, a, c, g, t, n}, 0, 2, []Base{T, A, T, A}, []Base{T, A, T, A, G, T, N, a, c, g, t, n}},
	{[]Base{A, C, G, T, N, a, c, g, t, n}, 8, 10, []Base{T, A, T, A}, []Base{A, C, G, T, N, a, c, g, T, A, T, A}},
	{[]Base{A, C, G, T, N, a, c, g, t, n}, 0, 10, []Base{A, C, G, T}, []Base{A, C, G, T}},
	{[]Base{A, C, G, T, N, a, c, g, t, n}, 3, 4, []Base{A, C, G, T}, []Base{A, C, G, A, C, G, T, N, a, c, g, t, n}},
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

func TestDel(t *testing.T) {
	for _, test := range deletionTests {
		actual := Delete(test.input, test.delStart, test.delEnd)
		if CompareSeqsCaseSensitive(actual, test.expected) != 0 {
			t.Errorf("Deleting positions %d to %d of %v gave %v when %v was expected.", test.delStart, test.delEnd, test.input, actual, test.expected)
		}
	}
}

func TestInsert(t *testing.T) {
	for _, test := range insertionTests {
		actual := Insert(test.seq, test.pos, test.insSeq)
		if CompareSeqsCaseSensitive(actual, test.expected) != 0 {
			t.Errorf("Inserting %v into %v at position %d gave %v when %v was expected.", test.insSeq, test.seq, test.pos, actual, test.expected)
		}
	}
}

func TestReplace(t *testing.T) {
	for _, test := range replaceTests {
		actual := Replace(test.seq, test.start, test.end, test.insSeq)
		if CompareSeqsCaseSensitive(actual, test.expected) != 0 {
			t.Errorf("Replacing %d to %d of %v with %v gave %v when %v was expected.", test.start, test.end, test.seq, test.insSeq, actual, test.expected)
		}
	}
}
