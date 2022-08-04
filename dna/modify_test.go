package dna

import "testing"

var revTests = []struct {
	input []Base
	expected []Base
}{
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, []Base{LowerN, LowerT, LowerG, LowerC, LowerA, N, T, G, C, A}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, []Base{N, T, G, C, A, N, T, G, C, A}},
}

var revCompTests = []struct {
	input    []Base // input
	expected []Base // after applying reverseComplement
}{
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, []Base{LowerN, LowerA, LowerC, LowerG, LowerT, N, A, C, G, T}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, []Base{N, A, C, G, T, N, A, C, G, T}},
	{[]Base{LowerA, LowerC, LowerG, LowerT, LowerN, A, C, G, T, N}, []Base{N, A, C, G, T, LowerN, LowerA, LowerC, LowerG, LowerT}},
	{[]Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}, []Base{LowerN, LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT}},
	{[]Base{A, A, A, G, G, T, T, T, N, LowerN, LowerC, LowerG, LowerT, LowerN}, []Base{LowerN, LowerA, LowerC, LowerG, LowerN, N, A, A, A, C, C, T, T, T}},
}

var deletionTests = []struct {
	input    []Base // input
	delStart int
	delEnd   int
	expected []Base // after applying deletion
}{
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, 0, 2, []Base{G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}},
	{[]Base{A, C, G, T, N, A, C, G, T, N}, 7, 10, []Base{A, C, G, T, N, A, C}},
	{[]Base{LowerA, LowerC, LowerG, LowerT, LowerN, A, C, G, T, N}, 3, 4, []Base{LowerA, LowerC, LowerG, LowerN, A, C, G, T, N}},
	{[]Base{LowerA, LowerC, LowerG, LowerT, LowerN, LowerA, LowerC, LowerG, LowerT, LowerN}, 4, 8, []Base{LowerA, LowerC, LowerG, LowerT, LowerT, LowerN}},
}

var insertionTests = []struct {
	seq      []Base // input
	pos      int
	insSeq   []Base
	expected []Base // after applying deletion
}{
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, 0, []Base{T, A, T, A}, []Base{T, A, T, A, A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}},
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, 10, []Base{T, A, T, A}, []Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN, T, A, T, A}},
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, 3, []Base{A, C, G, T}, []Base{A, C, G, A, C, G, T, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}},
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, 4, []Base{A, C, G, T}, []Base{A, C, G, T, A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}},
}

var replaceTests = []struct {
	seq      []Base // input
	start    int
	end      int
	insSeq   []Base
	expected []Base // after applying deletion
}{
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, 0, 2, []Base{T, A, T, A}, []Base{T, A, T, A, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}},
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, 8, 10, []Base{T, A, T, A}, []Base{A, C, G, T, N, LowerA, LowerC, LowerG, T, A, T, A}},
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, 0, 10, []Base{A, C, G, T}, []Base{A, C, G, T}},
	{[]Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}, 3, 4, []Base{A, C, G, T}, []Base{A, C, G, A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}},
}

func TestRev(t *testing.T) {
	var actual []Base
	for _, test := range revTests {
		actual = make([]Base, len(test.input))
		copy(actual, test.input)
		Reverse(actual)
		if CompareSeqsCaseSensitive(actual, test.expected) != 0 {
			t.Errorf("Reversing %v gave %v when %v was expected.", test.input, actual, test.expected)
		}
	}
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
	var actual []Base
	for _, test := range deletionTests {
		actual = Delete(test.input, test.delStart, test.delEnd)
		if CompareSeqsCaseSensitive(actual, test.expected) != 0 {
			t.Errorf("Deleting positions %d to %d of %v gave %v when %v was expected.", test.delStart, test.delEnd, test.input, actual, test.expected)
		}
	}
}

func TestInsert(t *testing.T) {
	var actual []Base
	for _, test := range insertionTests {
		actual = Insert(test.seq, test.pos, test.insSeq)
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

func TestComplement(t *testing.T) {
	expected := []Base{A, C, G, T, LowerA, LowerC, LowerG, LowerT, N, LowerN, Gap, Dot, Nil}
	input := []Base{T, G, C, A, LowerT, LowerG, LowerC, LowerA, N, LowerN, Gap, Dot, Nil}
	actual := []Base{T, G, C, A, LowerT, LowerG, LowerC, LowerA, N, LowerN, Gap, Dot, Nil}
	Complement(actual)
	if CompareSeqsCaseSensitive(actual, expected) != 0 {
		t.Errorf("Complementing %v gave %v when %v was expected.", input, actual, expected)
	}
}

func TestRemoveGaps(t *testing.T) {
	expected := []Base{A, C, G, T, LowerA, LowerC, LowerG, LowerT, N, LowerN, Dot, Nil}
	input := []Base{A, C, G, T, LowerA, LowerC, LowerG, LowerT, N, LowerN, Gap, Dot, Nil}
	actual := RemoveGaps(input)
	if CompareSeqsCaseSensitive(actual, expected) != 0 {
		t.Errorf("RemoveGaps on %v gave %v when %v was expected", input, actual, expected)
	}
}
