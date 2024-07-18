package cigar

import (
	"testing"
)

func TestAddCigar(t *testing.T) {
	tests := []struct {
		alpha    []Cigar
		beta     Cigar
		expected []Cigar
	}{
		{[]Cigar{{RunLength: 2, Op: Deletion}, {RunLength: 30, Op: Match}}, Cigar{RunLength: 20, Op: Match}, []Cigar{{RunLength: 2, Op: Deletion}, {RunLength: 50, Op: Match}}},
		{[]Cigar{{RunLength: 1, Op: Deletion}, {RunLength: 10, Op: Match}}, Cigar{RunLength: 2, Op: SoftClip}, []Cigar{{RunLength: 1, Op: Deletion}, {RunLength: 10, Op: Match}, {RunLength: 2, Op: SoftClip}}},
	}
	for _, test := range tests {
		result := AddCigar(test.alpha, test.beta)
		if !AllEqual(test.expected, result) {
			t.Errorf("Error: Incorrect AddCigar() %s != %s\n", ToString(result), ToString(test.expected))
		}
	}
}

func TestCatCigar(t *testing.T) {
	tests := []struct {
		alpha    []Cigar
		beta     []Cigar
		expected []Cigar
	}{
		{[]Cigar{}, []Cigar{{5, Match}, {3, Deletion}}, []Cigar{{5, Match}, {3, Deletion}}},                                         // Empty alpha
		{[]Cigar{{2, Insertion}}, []Cigar{}, []Cigar{{2, Insertion}}},                                                               // Empty beta
		{[]Cigar{{8, Match}}, []Cigar{{2, Match}}, []Cigar{{10, Match}}},                                                            // Merge at start
		{[]Cigar{{5, Deletion}}, []Cigar{{3, Insertion}, {4, 'X'}}, []Cigar{{5, Deletion}, {3, Insertion}, {4, 'X'}}},               // No merge
		{[]Cigar{{2, SoftClip}, {1, Match}}, []Cigar{{6, Match}, {1, Deletion}}, []Cigar{{2, SoftClip}, {7, Match}, {1, Deletion}}}, // Merge and append
	}

	for _, test := range tests {
		result := CatCigar(test.alpha, test.beta)
		if !AllEqual(test.expected, result) {
			t.Errorf("Error: Incorrect CatCigar() result. %s != %s", ToString(test.expected), ToString(result))
		}
	}
}

func TestReverseCigar(t *testing.T) {
	input := []Cigar{{7, Match}, {3, Deletion}, {3, Match}, {5, SoftClip}}
	expected := []Cigar{{5, SoftClip}, {3, Match}, {3, Deletion}, {7, Match}}

	// Make a copy for reversing to avoid modifying the original
	reversed := make([]Cigar, len(input))
	copy(reversed, input)
	ReverseCigar(reversed)

	if !AllEqual(reversed, expected) {
		t.Errorf("Error: Incorrect ReverseCigar() result. %s != %s", ToString(reversed), ToString(expected))
	}
}

func TestSoftClipBases(t *testing.T) {
	tests := []struct {
		front        int
		lengthOfRead int
		cig          []Cigar
		expected     []Cigar
	}{
		{0, 10, []Cigar{{10, Match}}, []Cigar{{10, Match}}},                                                           // No clipping
		{2, 10, []Cigar{{5, Match}, {5, Deletion}}, []Cigar{{2, SoftClip}, {5, Match}, {5, Deletion}, {3, SoftClip}}}, // Front clipping
		{0, 15, []Cigar{{5, Match}, {5, Deletion}}, []Cigar{{5, Match}, {5, Deletion}, {10, SoftClip}}},               // End clipping
		{2, 12, []Cigar{{5, Match}, {5, Deletion}}, []Cigar{{2, SoftClip}, {5, Match}, {5, Deletion}, {5, SoftClip}}}, // Both clippings
	}

	for _, test := range tests {
		result := SoftClipBases(test.front, test.lengthOfRead, test.cig)
		if !AllEqual(test.expected, result) {
			t.Errorf("Error: SoftClipBases(%d, %d, %s) = %s != %s", test.front, test.lengthOfRead, ToString(test.cig), ToString(result), ToString(test.expected))
		}
	}
}
