package dna

import (
	"testing"
)

var compareTests = []struct {
	a                               string // first input
	b                               string // second input
	expectedCaseSensitive           int
	expectedIgnoreCase              int
	expectedIgnoreCaseAndGaps       int
	expectedCaseSensitiveIgnoreGaps int
}{
	{"ACGTACGT", "ACGTacgt", -1, 0, 0, -1},
	{"ACGTACGTT", "ACGTacgt", -1, 1, 1, -1},
	{"ACGTACGT", "ACGTacgtT", -1, -1, -1, -1},
	{"ACGTA---CGT", "ACGTacgt", -1, 1, 0, -1},
	{"ACGTACGT", "ACGTac---gt", -1, -1, 0, -1},
	{"A-CGTACGT", "A---CGTacgt", -1, -1, 0, -1},
	{"ACGTACGT--", "ACGTacgt", -1, 1, 0, -1},
	{"GAGTGTGGATATACTTGCCTGTTCTGGGGGTGTACACGTG--TGTGTGCACACAGGTGACCCTGTACAGTGATTGCATGCGTGCACCAGGGAGTGTGGATATACTTGCCTGTTCTGGGGGTGTACA--", "GAGTGTGGATATACTTGCCTGTTCTGGGGGTGTACACGTGTGTGTGCACACAGGTGACCCTGTACAGTGATTGCATGCGTGCACCAGGGAGTGTGGATATACTTGCCTGTTCTGGGGGTGTACA", 1, 1, 0, 0},
}

func TestCompare(t *testing.T) {
	for _, test := range compareTests {
		a := StringToBases(test.a)
		b := StringToBases(test.b)
		answer := CompareSeqsCaseSensitive(a, b)
		if answer != test.expectedCaseSensitive {
			t.Errorf("CompareSeqsCaseSensitive: expected %d, got %d when testing %v and %v\n", test.expectedCaseSensitive, answer, a, b)
		}
		answer = CompareSeqsIgnoreCase(a, b)
		if answer != test.expectedIgnoreCase {
			t.Errorf("CompareSeqsIgnoreCase: expected %d, got %d when testing %v and %v\n", test.expectedIgnoreCase, answer, a, b)
		}
		answer = CompareSeqsIgnoreCaseAndGaps(a, b)
		if answer != test.expectedIgnoreCaseAndGaps {
			t.Errorf("CompareSeqsIgnoreCaseAndGaps: expected %d, got %d when testing %v and %v\n", test.expectedIgnoreCaseAndGaps, answer, a, b)
		}
		answer = CompareSeqsCaseSensitiveIgnoreGaps(a, b)
		if answer != test.expectedCaseSensitiveIgnoreGaps {
			t.Errorf("CompareSeqsCaseSensitiveIgnoreGaps: expected %d, got %d when testing %v and %v\n", test.expectedCaseSensitiveIgnoreGaps, answer, a, b)
		}
	}
}

//TODO Compare2DSeqsTests
