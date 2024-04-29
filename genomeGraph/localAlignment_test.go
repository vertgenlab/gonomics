package genomeGraph

import (
	"log"
	"testing"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
)

// func RightLocal(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]rune) (int64, []*cigar.Cigar, int64, int64, int64, int64)
func swMatrixSetup(size int64) ([][]int64, [][]byte) {
	m := make([][]int64, size)
	trace := make([][]byte, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		trace[idx] = make([]byte, size)
	}
	return m, trace
}

func TestRightLocal(t *testing.T) {
	var seqOneA = dna.StringToBases("TTTTTTTTTTTTTTTTAGC")
	var seqOneB = dna.StringToBases("ATTTTTTTTTTTTTTTTAGC")
	m, trace := swMatrixSetup(10000)

	score, alignmentPath, refStart, refEnd, queryStart, queryEnd := RightLocal(seqOneA, seqOneB, align.HumanChimpTwoScoreMatrix, -600, m, trace)

	log.Printf("score=%d, alignment=%s, refStart=%d, refEnd=%d, queryStart=%d, queryEnd=%d\n", score, cigar.ByteCigarToString(alignmentPath), refStart, refEnd, queryStart, queryEnd)
}

func BenchmarkLocalLeftDna(b *testing.B) {
	var seqOneA = dna.StringToBases("TAGGGGGTGGGGGGGGT")
	var seqOneB = dna.StringToBases("CAGGGGGTGGGGGGGG")
	m, trace := swMatrixSetup(10000)
	b.ResetTimer() // Resets the timer to exclude setup time from benchmark measurements
	b.ReportAllocs()
	for n := 0; n < b.N; n++ {
		LeftLocal(seqOneA, seqOneB, align.HumanChimpTwoScoreMatrix, -600, m, trace)
	}
	b.StopTimer()
}

// BenchmarkTwoBitLocalLeftAlign benchmarks the TwoBitLocalLeftAlign function.
func BenchmarkTwoBitLocalLeftAlign(b *testing.B) {
	var seqOneA = dnaTwoBit.NewTwoBit(dna.StringToBases("TAGGGGGTGGGGGGGGT"))
	var seqOneB = dnaTwoBit.NewTwoBit(dna.StringToBases("CAGGGGGTGGGGGGGG"))
	m, trace := swMatrixSetup(10000)
	b.ReportAllocs()
	b.ResetTimer() // Resets the timer to exclude setup time from benchmark measurements
	for n := 0; n < b.N; n++ {
		TwoBitLocalLeftAlign(seqOneA, seqOneB, align.HumanChimpTwoScoreMatrix, -600, m, trace)
	}
	b.StopTimer()

}
