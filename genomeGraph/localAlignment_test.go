package genomeGraph

import (
	"testing"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

var tests = []struct {
	Alpha                                                  []dna.Base
	Beta                                                   []dna.Base
	ExpectedScore                                          int64
	ExpectedCigar                                          []cigar.Cigar
	ExpectedMinI, ExpectedMaxI, ExpectedMinJ, ExpectedMaxJ int
}{
	{
		Alpha:         dna.StringToBases("TTTTTTTTTTTTTTTTAGC"),
		Beta:          dna.StringToBases("TTTTTTAGC"),
		ExpectedScore: int64(830),
		ExpectedCigar: []cigar.Cigar{{RunLength: 9, Op: '='}},
		ExpectedMinI:  10,
		ExpectedMaxI:  19,
		ExpectedMinJ:  0,
		ExpectedMaxJ:  9,
	},
	{
		Alpha:         dna.StringToBases("GACCCTGACCTTACTAGTTTACAATCACACGATCCTGACCTTACTAGTTT"),
		Beta:          dna.StringToBases("GACCCTGACCTCACACGATCCTGACCTTACTAGTTT"),
		ExpectedScore: int64(2450),
		ExpectedCigar: []cigar.Cigar{{RunLength: 26, Op: '='}},
		ExpectedMinI:  24,
		ExpectedMaxI:  50,
		ExpectedMinJ:  10,
		ExpectedMaxJ:  36,
	},
}

var alphaRight = dna.StringToBases("TTTTTTTTTTTTTTTTAGC")
var betaRight = dna.StringToBases("TTTTTTTTTAGC")
var alphaLeft = dna.StringToBases("TAGGGGGTGGGGGGGGT")
var betaLeft = dna.StringToBases("GGGGGGGT")

func TestRightLocal(t *testing.T) {
	config := DefaultAlignment()
	pool := NewMatrixPool(50)

	score, alignmentPath, refStart, refEnd, queryStart, queryEnd := RightLocal(alphaRight, betaRight, config, pool)
	t.Logf("RightLocal:\n\n%s\nscore=%d, alignment=%s, refStart=%d, refEnd=%d, queryStart=%d, queryEnd=%d\n", View(alphaRight, betaRight, alignmentPath, refStart, refEnd, queryStart, queryEnd), score, cigar.ToString(alignmentPath), refStart, refEnd, queryStart, queryEnd)

	if score != 810 {
		t.Errorf("Score mismatch: expected %d, got %d", score, 810)
	}
	if !cigar.AllEqual(alignmentPath, []cigar.Cigar{{RunLength: 9, Op: '='}}) {
		t.Errorf("Alignment mismatch:\nExpected: %s\nGot: %s", cigar.ToString(alignmentPath), cigar.ToString([]cigar.Cigar{{RunLength: 9, Op: '='}}))
	}
	if refStart != 0 {
		t.Errorf("refStart mismatch: expected %d, got %d", 0, refStart)
	}
	if refEnd != 9 {
		t.Errorf("refEnd mismatch: expected %d, got %d", 9, refEnd)
	}
	if queryStart != 0 {
		t.Errorf("queryStart mismatch: expected %d, got %d", 0, queryStart)
	}
	if queryEnd != 9 {
		t.Errorf("queryEnd mismatch: expected %d, got %d", 9, queryEnd)
	}
}

func TestLeftLocal(t *testing.T) {
	config := DefaultAlignment()
	pool := NewMatrixPool(50)
	score, alignmentPath, refStart, refEnd, queryStart, queryEnd := LeftLocal(alphaLeft, betaLeft, config, pool)
	t.Logf("LeftLocal:\n\n%s\nscore=%d, alignment=%s, refStart=%d, refEnd=%d, queryStart=%d, queryEnd=%d\n", View(alphaLeft, betaLeft, alignmentPath, refStart, refEnd, queryStart, queryEnd), score, cigar.ToString(alignmentPath), refStart, refEnd, queryStart, queryEnd)
	if score != 790 {
		t.Errorf("Score mismatch: expected %d, got %d", score, 790)
	}
	if !cigar.AllEqual(alignmentPath, []cigar.Cigar{{RunLength: 8, Op: cigar.Equal}}) {
		t.Errorf("Alignment mismatch:\nExpected: %s\nGot: %s", cigar.ToString(alignmentPath), cigar.ToString([]cigar.Cigar{{RunLength: 8, Op: '='}}))
	}
	if refStart != 9 {
		t.Errorf("refStart mismatch: expected %d, got %d", 9, refStart)
	}
	if refEnd != 17 {
		t.Errorf("refEnd mismatch: expected %d, got %d", 17, refEnd)
	}
	if queryStart != 0 {
		t.Errorf("queryStart mismatch: expected %d, got %d", 0, queryStart)
	}
	if queryEnd != 8 {
		t.Errorf("queryEnd mismatch: expected %d, got %d", 8, queryEnd)
	}
}

func TestSmithWaterman(t *testing.T) {
	config := DefaultAlignment()
	dp := NewMatrixPool(20)
	for _, test := range tests {
		score, route, minI, maxI, minJ, maxJ := SmithWaterman(test.Alpha, test.Beta, config, dp)
		name := View(test.Alpha, test.Beta, route, minI, maxI, minJ, maxJ)
		if score != test.ExpectedScore {
			t.Errorf("Test Case:\n\n%s\nExpected score %d, got %d\n", name, test.ExpectedScore, score)
		}
		if !cigar.AllEqual(route, test.ExpectedCigar) {
			t.Errorf("Test Case:\n\n%s\nExpected cigar %s != %s\n", name, cigar.ToString(test.ExpectedCigar), cigar.ToString(route))
		}
		if minI != test.ExpectedMinI || maxI != test.ExpectedMaxI || minJ != test.ExpectedMinJ || maxJ != test.ExpectedMaxJ {
			t.Errorf("Test Case:\n\n%s\nExpected indices mismatch: got minI=%d, maxI=%d, minJ=%d, maxJ=%d\n", name, minI, maxI, minJ, maxJ)
		}
	}
}

func BenchmarkSmithWatermanDP(b *testing.B) {
	config := DefaultAlignment()
	pool := NewMatrixPool(500)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for _, test := range tests {
			score, route, minI, maxI, minJ, maxJ := SmithWaterman(test.Alpha, test.Beta, config, pool)
			if score != test.ExpectedScore {
				b.Errorf("Expected score %d, got %d\n", test.ExpectedScore, score)
			}
			if !cigar.AllEqual(route, test.ExpectedCigar) {
				b.Errorf("Expected cigar %s != %s\n", cigar.ToString(test.ExpectedCigar), cigar.ToString(route))
			}
			if minI != test.ExpectedMinI || maxI != test.ExpectedMaxI || minJ != test.ExpectedMinJ || maxJ != test.ExpectedMaxJ {
				b.Errorf("Expected indices mismatch: got minI=%d, maxI=%d, minJ=%d, maxJ=%d\n", minI, maxI, minJ, maxJ)
			}
		}
	}
	b.ReportAllocs()
}

func BenchmarkLeftLocalAlignment(b *testing.B) {
	config := DefaultAlignment()
	pool := NewMatrixPool(500)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		score, route, minI, maxI, minJ, maxJ := LeftLocal(alphaLeft, betaLeft, config, pool)
		b.Logf("Test case: Left local alignment\n%s\n", View(alphaLeft, betaLeft, route, minI, maxI, minJ, maxJ))
		if score != 790 {
			b.Errorf("Expected score %d, got %d\n", 790, score)
		}
	}
	b.ReportAllocs()
}

func BenchmarkRightLocalAlignment(b *testing.B) {
	config := DefaultAlignment()
	pool := NewMatrixPool(500)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		score, route, minI, maxI, minJ, maxJ := RightLocal(alphaRight, betaRight, config, pool)
		b.Logf("Test case: Right local alignment\n%s\n", View(alphaRight, betaRight, route, minI, maxI, minJ, maxJ))
		if score != 810 {
			b.Errorf("Score mismatch: expected %d, got %d", score, 810)
		}
	}
	b.ReportAllocs()
}
