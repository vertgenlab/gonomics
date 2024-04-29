package genomeGraph

import (
	"log"
	"testing"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
)

func TestGetLeftTwoBitBases(t *testing.T) {
	n1 := &Node{
		Seq: dna.StringToBases("ACGTA"),
	}
	n1.SeqTwoBit = dnaTwoBit.NewTwoBit(n1.Seq)

	n2 := &Node{
		Seq: dna.StringToBases("ACGATCGAT"),
	}
	n2.SeqTwoBit = dnaTwoBit.NewTwoBit(n2.Seq)

	AddEdge(n1, n2, 1)

	right := dna.StringToBases("AC")
	rightExpected := dna.StringToBases("ACGTA")
	rightTargetBases := RightGetTwoBit(n1, 10, 2, dnaTwoBit.NewTwoBit(right), dnaTwoBit.NewTwoBit(dna.StringToBases("")))

	if dnaTwoBit.ToString(rightTargetBases) != dna.BasesToString(rightExpected) {
		t.Errorf("getRightBases() = %v, want %v", dnaTwoBit.ToString(rightTargetBases), dna.BasesToString(rightExpected))
	}
	leftSeq := dna.StringToBases("A")
	leftExpected := dna.StringToBases("CGATCGA")
	leftTargetBases := LeftGetTwoBit(n2, 7, 7, dnaTwoBit.NewTwoBit(leftSeq), dnaTwoBit.NewTwoBit(dna.StringToBases("")))

	if dna.BasesToString(leftExpected) != dnaTwoBit.ToString(leftTargetBases) {
		t.Errorf("getLeftTargetBases() = %v, want %v", dna.BasesToString(leftExpected), dnaTwoBit.ToString(leftTargetBases))
	}
}

func TestTwoBitLeftLocal(t *testing.T) {
	var seqOneA = dna.StringToBases("TAGGGGGTGGGGGGGGT")
	var seqOneB = dna.StringToBases("CAGGGGGTGGGGGGGG")
	m, trace := swMatrixSetup(10000)

	score, alignmentPath, refStart, refEnd, queryStart, queryEnd := TwoBitLocalLeftAlign(dnaTwoBit.NewTwoBit(seqOneA), dnaTwoBit.NewTwoBit(seqOneB), align.HumanChimpTwoScoreMatrix, -600, m, trace)
	expextedScore := int64(880)
	log.Printf("score=%d, alignment=%s, refStart=%d, refEnd=%d, queryStart=%d, queryEnd=%d\n", score, cigar.ByteCigarToString(alignmentPath), refStart, refEnd, queryStart, queryEnd)
	if score != expextedScore {
		t.Errorf("Error: Left direction expected score %d, got %d", expextedScore, score)
	}
}

func TestLeftTwoBitLocal(t *testing.T) {
	var seqOneA = dnaTwoBit.NewTwoBit(dna.StringToBases("TAGGGGGTGGGGGGGGT"))
	var seqOneB = dnaTwoBit.NewTwoBit(dna.StringToBases("CAGGGGGTGGGGGGGG"))
	m, trace := swMatrixSetup(10000)
	score, alignmentPath, refStart, refEnd, queryStart, queryEnd := TwoBitLocalLeftAlign(seqOneA, seqOneB, align.HumanChimpTwoScoreMatrix, -600, m, trace)
	log.Printf("score=%d, alignment=%s, refStart=%d, refEnd=%d, queryStart=%d, queryEnd=%d\n", score, cigar.ByteCigarToString(alignmentPath), refStart, refEnd, queryStart, queryEnd)
}

// if rightScore != expectedRightScore {
// 		t.Errorf("Error: Right direction expected score %d, got %d", expectedRightScore, rightScore)
// 	}

// 	if cigar.ByteCigarToString(rightAlign) != cigar.ByteCigarToString(expectedRightAlignment) {
// 		t.Errorf("Error: Right direction expected alignment %+v, got %+v", expectedRightAlignment, rightAlign)
// 	}
// 	if cigar.ByteCigarToString(leftAlign) != cigar.ByteCigarToString(expectedLeftAlignment) {
// 		t.Errorf("Error: Left direction expected alignment %+v, got %+v", expectedLeftAlignment, leftAlign)
// 	}
// 	if expectedRightTargetEnd != rightTargetEnd || rightQueryEnd != expectedRightQueryEnd {
// 		t.Errorf("Error: Right direction: expected target and query start (%d, %d), got (%d, %d)", expectedRightTargetEnd, rightQueryEnd, rightTargetEnd, rightQueryEnd)
// 	}
// 	if expectedLeftTargetStart != leftTargetStart || leftQueryStart != expectedLeftQueryStart {
// 		t.Errorf("Error: Left direction expected target and query start (%d, %d), got (%d, %d)", expectedLeftTargetStart, leftQueryStart, leftTargetStart, leftQueryStart)
// 	}
