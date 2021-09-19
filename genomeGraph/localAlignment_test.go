package genomeGraph

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"testing"
)

//func RightLocal(alpha []dna.Base, beta []dna.Base, scores [][]int, gapPen int, m [][]int, trace [][]rune) (int, []*cigar.Cigar, int, int, int, int)

func TestRightLocal(t *testing.T) {
	var seqOneA = dna.StringToBases("TTTTTTTTTTTTTTTTAGC")
	var seqOneB = dna.StringToBases("ATTTTTTTTTTTTTTTTAGC")
	m, trace := swMatrixSetup(10000)

	score, alignmentPath, refStart, refEnd, queryStart, queryEnd := RightLocal(seqOneA, seqOneB, HumanChimpTwoScoreMatrix, -600, m, trace)

	log.Printf("score=%d, alignment=%s, refStart=%d, refEnd=%d, queryStart=%d, queryEnd=%d\n", score, cigar.ToString(alignmentPath), refStart, refEnd, queryStart, queryEnd)
}

func TestLeftLocal(t *testing.T) {
	var seqOneA = dna.StringToBases("TAGGGGGTGGGGGGGGT")
	var seqOneB = dna.StringToBases("CAGGGGGTGGGGGGGG")
	m, trace := swMatrixSetup(10000)

	score, alignmentPath, refStart, refEnd, queryStart, queryEnd := LeftLocal(seqOneA, seqOneB, HumanChimpTwoScoreMatrix, -600, m, trace)

	log.Printf("score=%d, alignment=%s, refStart=%d, refEnd=%d, queryStart=%d, queryEnd=%d\n", score, cigar.ToString(alignmentPath), refStart, refEnd, queryStart, queryEnd)
}
