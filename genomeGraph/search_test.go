package genomeGraph

import (
	"math/rand"
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
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

func leftBases(n *Node, extension int, index int, ans []dna.Base) []dna.Base {
	seqLen := len(ans)
	var basesToTake int = (numbers.Min(seqLen+index, extension) - seqLen)

	ans = append(ans, make([]dna.Base, basesToTake)...)
	copy(ans[seqLen:], n.Seq[index-basesToTake:index])
	return ans
}

func BenchmarkNewGetLeftBases(b *testing.B) {
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
		leftBases(n, 200, 500, ans) // Example parameters
	}
	b.ReportAllocs()
}

// func BenchmarkGetRightBases(b *testing.B) {
//      n := &Node{Seq: make([]dna.Base, 1000)} // Example 1000-base sequence
//     for i := range n.Seq {
//         n.Seq[i] = dna.Base(rand.Intn(4)) // Random bases
//     }

//     seq := make([]dna.Base, 100)        // Example 100-base query
//     ans := make([]dna.Base, 0, 200)     // Preallocate to avoid reallocations

//     b.ResetTimer()

//     for i := 0; i < b.N; i++ {
//         getRightBases(n, 200, 500, append(ans, seq...))
//     }
// }
