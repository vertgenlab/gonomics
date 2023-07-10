package cigar

import (
	"strings"
	"testing"
)

func TestFromString(t *testing.T) {
	var cigarsString = "35M2I16D"
	var c1 Cigar = Cigar{RunLength: 35, Op: 'M'}
	var c2 Cigar = Cigar{RunLength: 2, Op: 'I'}
	var c3 Cigar = Cigar{RunLength: 16, Op: 'D'}
	var cigars []Cigar = []Cigar{c1, c2, c3}

	cigarCheck := FromString(cigarsString)

	for i := 0; i < len(cigarCheck); i++ {
		if !isEqual(cigars[i], cigarCheck[i]) {
			t.Errorf("Error with FromString")
		}
	}
}

func TestToString(t *testing.T) {
	var cigarsString = "35M2I16D"
	var c1 Cigar = Cigar{RunLength: 35, Op: 'M'}
	var c2 Cigar = Cigar{RunLength: 2, Op: 'I'}
	var c3 Cigar = Cigar{RunLength: 16, Op: 'D'}
	var cigars []Cigar = []Cigar{c1, c2, c3}

	cigarCheck := ToString(cigars)

	if strings.Compare(cigarCheck, cigarsString) != 0 {
		t.Errorf("Error with ToString")
	}
}

func isEqual(a Cigar, b Cigar) bool {
	return (a.RunLength == b.RunLength && a.Op == b.Op)
}
