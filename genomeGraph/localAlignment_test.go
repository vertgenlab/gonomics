package genomeGraph

import (
	"testing"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

func TestRightLocal(t *testing.T) {
	var seqOneA = dna.StringToBases("TTTTTTTTTTTTTTTTAGC")
	var seqOneB = dna.StringToBases("ATTTTTTTTTTTTTTTTAGC")
	m, trace := swMatrixSetup(10000)

	score, alignmentPath, refStart, refEnd, queryStart, queryEnd := RightLocal(seqOneA, seqOneB, align.HumanChimpTwoScoreMatrix, -600, m, trace)

	t.Logf("score=%d, alignment=%s, refStart=%d, refEnd=%d, queryStart=%d, queryEnd=%d\n", score, cigar.ToString(alignmentPath), refStart, refEnd, queryStart, queryEnd)
}

func TestLeftLocal(t *testing.T) {
	var seqOneA = dna.StringToBases("TAGGGGGTGGGGGGGGT")
	var seqOneB = dna.StringToBases("CAGGGGGTGGGGGGGG")
	m, trace := swMatrixSetup(10000)

	score, alignmentPath, refStart, refEnd, queryStart, queryEnd := LeftLocal(seqOneA, seqOneB, align.HumanChimpTwoScoreMatrix, -600, m, trace)

	t.Logf("score=%d, alignment=%s, refStart=%d, refEnd=%d, queryStart=%d, queryEnd=%d\n", score, cigar.ToString(alignmentPath), refStart, refEnd, queryStart, queryEnd)
}
