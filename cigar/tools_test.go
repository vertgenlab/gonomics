package cigar

import (
	"testing"
)

func TestAppend(t *testing.T) {
	tests := []struct {
		alpha    []Cigar
		beta     Cigar
		expected []Cigar
	}{
		{[]Cigar{{2, Deletion}, {30, Match}}, Cigar{20, Match}, []Cigar{{2, Deletion}, {50, Match}}},
		{[]Cigar{{1, Deletion}, {10, Match}}, Cigar{2, SoftClip}, []Cigar{{1, Deletion}, {10, Match}, {2, SoftClip}}},
	}
	for _, test := range tests {
		result := Append(test.alpha, test.beta)
		if !AllEqual(test.expected, result) {
			t.Errorf("Error: Incorrect AddCigar() %s != %s\n", ToString(result), ToString(test.expected))
		}
	}
}

func TestConcat(t *testing.T) {
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
		result := Concat(test.alpha, test.beta)
		if !AllEqual(test.expected, result) {
			t.Errorf("Error: Incorrect CatCigar() result. %s != %s", ToString(test.expected), ToString(result))
		}
	}
}

func TestAppendSoftClips(t *testing.T) {
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
		result := AppendSoftClips(test.front, test.lengthOfRead, test.cig)
		if !AllEqual(test.expected, result) {
			t.Errorf("Error: AppendSoftClips(%d, %d, %s) = %s != %s", test.front, test.lengthOfRead, ToString(test.cig), ToString(result), ToString(test.expected))
		}
	}
}

func TestReverseCigar(t *testing.T) {
	tests := []struct {
		input    []Cigar
		expected []Cigar
	}{
		{[]Cigar{{7, Match}, {3, Deletion}, {3, Match}, {5, SoftClip}}, []Cigar{{5, SoftClip}, {3, Match}, {3, Deletion}, {7, Match}}},
		{[]Cigar{{7, Match}, {3, Deletion}, {3, Match}, {2, Mismatch}, {5, SoftClip}}, []Cigar{{5, SoftClip}, {2, Mismatch}, {3, Match}, {3, Deletion}, {7, Match}}}, // Odd number
	}
	for _, test := range tests {
		// Make a copy for reversing to avoid modifying the original
		reversed := make([]Cigar, len(test.input))
		copy(reversed, test.input)
		ReverseCigar(reversed)
		if !AllEqual(test.expected, reversed) {
			t.Errorf("Error: Incorrect ReverseCigar() result. %s != %s", ToString(reversed), ToString(test.expected))
		}
	}
}

func TestToUint32(t *testing.T) {
	var input Cigar = Cigar{10, Match}
	var expected uint32 = 160

	result := ToUint32(input)
	if result != expected {
		t.Errorf("Error: Incorrect ToUint32() result. %d != %d", result, expected)
	}
}

func TestIsUnmapped(t *testing.T) {
	tests := []struct {
		cigar    []Cigar
		expected bool
	}{
		{[]Cigar{}, true},
		{[]Cigar{{1, Match}, {1, Equal}, {1, Mismatch}}, false},
		{make([]Cigar, 0), true},
		{nil, true}, // a slice that is nil will currently return true 0xff
		{FromString("*"), true},
		{FromString("0xff"), true},
		{FromString("150M"), false},
		{make([]Cigar, 1), false}, // TODO: Consider the result in which a slice is allocated memory, but contains empty values: i.e Cigar{0, nil}.
	}
	for _, test := range tests {
		result := IsUnmapped(test.cigar)
		if test.expected != result {
			t.Errorf("Error: Incorrect IsUnmapped() result. %s (%t) != %t", ToString(test.cigar), result, test.expected)
		}
	}
}
