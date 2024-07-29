package cigar

import (
	"testing"
)

func TestAllEqual(t *testing.T) {
	tests := []struct {
		alpha    []Cigar
		beta     []Cigar
		expected bool
	}{
		{[]Cigar{{RunLength: 10, Op: Match}}, []Cigar{{RunLength: 2, Op: Deletion}, {RunLength: 5, Op: Match}}, false},
		{[]Cigar{{RunLength: 150, Op: Match}}, []Cigar{{RunLength: 150, Op: Match}}, true},
		{[]Cigar{{RunLength: 1, Op: Deletion}}, []Cigar{{RunLength: 1, Op: Insertion}}, false},
	}
	for _, test := range tests {
		result := AllEqual(test.alpha, test.beta)
		if test.expected != result {
			t.Errorf("Error: Incorrect AllEqual() %v != %v\n", result, test.expected)
		}
	}
}
