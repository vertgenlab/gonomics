package genomeGraph

import (
	"math/rand"
	"reflect"
	"sync"
	"testing"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

func TestGetTargetBases(t *testing.T) {
	n1 := &Node{
		Seq: dna.StringToBases("ACGTA"),
	}
	n2 := &Node{
		Seq: dna.StringToBases("ACGATCGAT"),
	}
	AddEdge(n1, n2, 1)

	rightSeq := dna.StringToBases("AC")
	rightExpected := dna.StringToBases("ACGTA")
	rightTargetBases := getRightBases(n1, 10, 2, rightSeq)

	if dna.BasesToString(rightTargetBases) != dna.BasesToString(rightExpected) {
		t.Errorf("getRightBases() = %v, want %v", dna.BasesToString(rightTargetBases), dna.BasesToString(rightExpected))
	}
	leftSeq := dna.StringToBases("A")
	leftExpected := dna.StringToBases("ACGATCG")
	leftTargetBases := getLeftBases(n2, 7, 7, leftSeq)

	if dna.BasesToString(leftExpected) != dna.BasesToString(leftTargetBases) {
		t.Errorf("getLeftTargetBases() = %v, want %v", dna.BasesToString(leftTargetBases), dna.BasesToString(leftExpected))
	}

}
func TestGetLeftBases(t *testing.T) {
	n := &Node{Seq: dna.StringToBases("ACGTAGCACGTAGCACGTAGC")} // Longer sequence

	tests := []struct {
		name      string
		extension int
		refEnd    int // Changed from index to refEnd
		seq       []dna.Base
		ans       []dna.Base
		expected  []dna.Base
	}{
		{
			name:      "Basic extension",
			extension: 3,
			refEnd:    3, // Changed index to refEnd
			ans:       dna.StringToBases("GC"),
			expected:  dna.StringToBases("GCG"), // Expected output based on new refEnd
		},
		{
			name:      "Extension exceeds available bases",
			extension: 5,
			refEnd:    2, // Changed index to refEnd
			ans:       dna.StringToBases("A"),
			expected:  dna.StringToBases("AAC"), // Expected output based on new refEnd
		},
		{
			name:      "No extension needed",
			extension: 2,
			refEnd:    4, // Changed index to refEnd
			ans:       dna.StringToBases("GC"),
			expected:  dna.StringToBases("GC"),
		},
		{
			name:      "Empty seq and ans",
			extension: 4,
			refEnd:    3, // Changed index to refEnd
			seq:       dna.StringToBases(""),
			ans:       dna.StringToBases(""),
			expected:  dna.StringToBases("ACG"),
		},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			result := getLeftBases(n, test.extension, test.refEnd, test.ans)
			if dna.BasesToString(result) != dna.BasesToString(test.expected) {
				t.Errorf("getLeftBases(%d, %d, %q, %q) = %q; want %q",
					test.extension, test.refEnd, dna.BasesToString(test.seq), dna.BasesToString(test.seq),
					dna.BasesToString(result), dna.BasesToString(test.expected))
			}
		})
	}
}

func BenchmarkGetLeftBases(b *testing.B) {
	// Create test data
	n := &Node{Seq: make([]dna.Base, 1000)} // Example 1000-base sequence
	for i := range n.Seq {
		n.Seq[i] = dna.Base(rand.Intn(4)) // Random bases
	}

	seq := make([]dna.Base, 100)    // Example 100-base query
	ans := make([]dna.Base, 0, 200) // Preallocate to avoid reallocations
	ans = append(ans, seq...)
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		getLeftBases(n, 200, 500, ans) // Example parameters
	}
	b.ReportAllocs()
}

func TestRightAlignTraversal(t *testing.T) {
    // Test Data Setup
    config := DefaultAlignment()
	matrix := NewMatrixPool(defaultMatrixSize)
    pool := &sync.Pool{
        New: func() interface{} { return &dnaPool{} },
    }

    // Create sample nodes, sequences, and reads
    n1 := &Node{Seq: dna.StringToBases("ACGT")}
	n2 := &Node{Seq: dna.StringToBases("TTG")}
	AddEdge(n1, n2, 1)
    seq :=  dna.StringToBases("ACTG")
    read :=  dna.StringToBases("TG")

    // Test Cases
    testCases := []struct {
        name       string
        start      int
        extension  int
        expected   []cigar.Cigar
        expectScore int64
    }{
        {"Simple Path", 0, 2, []cigar.Cigar{{Op: 'M', RunLength: 2}}, 10}, // Adjust expected score as needed
        // ... add more test cases with different scenarios ...
    }

    for _, tc := range testCases {
        t.Run(tc.name, func(t *testing.T) {
            sk := scoreKeeper{}
            dynamicScore := dynamicScoreKeeper{}

            // Call the function under test
            result, score, _, _, _ := RightAlignTraversal(n1, seq, tc.start, []uint32{}, read, tc.extension, config, matrix, sk, dynamicScore, pool)

            // Assertions
            if !reflect.DeepEqual(result, tc.expected) {
                t.Errorf("Expected alignment %v, got %v", tc.expected, result)
            }
            if score != tc.expectScore {
                t.Errorf("Expected score %d, got %d", tc.expectScore, score)
            }

        })
    }
}
