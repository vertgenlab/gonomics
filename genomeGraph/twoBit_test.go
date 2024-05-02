package genomeGraph

import (
	"testing"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
)

func TestGetwoBitBases(t *testing.T) {
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

	rightDnaBases := getRightBases(n1, 10, 2, right, dna.StringToBases(""))

	if dnaTwoBit.ToString(rightTargetBases) != dna.BasesToString(rightDnaBases) {
		t.Errorf("RightGetTwoBit() = %v, want %v", dna.BasesToString(rightDnaBases), dnaTwoBit.ToString(rightTargetBases))
	}

	if dnaTwoBit.ToString(rightTargetBases) != dna.BasesToString(rightExpected) {
		t.Errorf("getRightBases() = %v, want %v", dnaTwoBit.ToString(rightTargetBases), dna.BasesToString(rightExpected))
	}
	leftSeq := dna.StringToBases("A")
	leftExpected := dna.StringToBases("CGATCGA")
	leftTargetBases := LeftGetTwoBit(n2, 7, 7, dnaTwoBit.NewTwoBit(leftSeq), dnaTwoBit.NewTwoBit(dna.StringToBases("")))

	leftDnaBases := getLeftBases(n2, 7, 7, leftSeq, dna.StringToBases(""))

	if dnaTwoBit.ToString(leftTargetBases) != dna.BasesToString(leftDnaBases) {
		t.Errorf("LeftGetTwoBit() = %v, want %v", dna.BasesToString(leftDnaBases), dnaTwoBit.ToString(leftTargetBases))
	}
	if dna.BasesToString(leftExpected) != dnaTwoBit.ToString(leftTargetBases) {
		t.Errorf("LeftGetTwoBit()= %v, want %v", dna.BasesToString(leftExpected), dnaTwoBit.ToString(leftTargetBases))
	}
}

func TestTwoBitLeftLocal(t *testing.T) {
	var seqOneA = dna.StringToBases("TAGGGGGTGGGGGGGGT")
	var seqOneB = dna.StringToBases("CAGGGGGTGGGGGGGG")
	m, trace := swMatrixSetup(10000)
	config := &GraphSettings{
		ScoreMatrix:    align.HumanChimpTwoScoreMatrix,
		TileSize:       32,
		StepSize:       32,
		GapPenalty:     -600,
		OpenGapPenalty: -150,
	}
	memoryPool := MatrixPoolMemory(defaultMatrixSize)

	scoreTwoBit, alignmentPathTwoBit, targetStartTwoBit, queryStartTwoBit := TwoBitLocalLeftAlign(dnaTwoBit.NewTwoBit(seqOneA), dnaTwoBit.NewTwoBit(seqOneB), config, memoryPool)
	scoreDnaBase, alignmentPathDnaBase, targetStartDnaBase, targetEndDnaBase, queryStartDnaBase, queryEndDnaBase := LeftLocal(seqOneA, seqOneB, align.HumanChimpTwoScoreMatrix, config.GapPenalty, m, trace)

	expextedScore := int64(880)
	expectedAlign := "15=1D"
	expectedTargetStart := 1
	expectedTargetEnd := 17
	expectedQueryStart := 1
	expectedQueryEnd := 16
	t.Logf("score=%d, alignment=%s, refStart=%d, queryStart=%d\n", scoreTwoBit, cigar.ByteCigarToString(alignmentPathTwoBit), targetStartTwoBit, queryStartTwoBit)

	if scoreTwoBit != expextedScore || scoreTwoBit != scoreDnaBase {
		t.Errorf("Error: Left direction expected score %d, got %d", scoreDnaBase, scoreTwoBit)
	}
	if expectedAlign != cigar.ByteCigarToString(alignmentPathDnaBase) || cigar.ByteCigarToString(alignmentPathTwoBit) != cigar.ByteCigarToString(alignmentPathDnaBase) {
		t.Errorf("Error: Left direction expected score %s, got %s", expectedAlign, cigar.ByteCigarToString(alignmentPathTwoBit))
	}

	if targetStartTwoBit != expectedTargetStart || targetStartTwoBit != targetStartDnaBase {
		t.Errorf("Error: Left target start position %d, != %d", expectedTargetStart, targetStartTwoBit)
	}

	if expectedTargetEnd != targetEndDnaBase {
		t.Errorf("Error: Left target end position %d, != %d", expectedTargetEnd, targetEndDnaBase)
	}

	if queryStartTwoBit != expectedQueryStart || queryStartTwoBit != queryStartDnaBase {
		t.Errorf("Error: Left query start position %d, != %d", expectedQueryStart, scoreTwoBit)
	}

	if expectedQueryEnd != queryEndDnaBase {
		t.Errorf("Error: Left query end position %d, != %d", expectedQueryEnd, queryEndDnaBase)
	}
}
func TestTwoBitRightLocal(t *testing.T) {
	var seqOneA = dna.StringToBases("TAGGGGGTGGGGGGGGT")
	var seqOneB = dna.StringToBases("CAGGGGGTGGGGGGGG")

	config := &GraphSettings{
		ScoreMatrix:    align.HumanChimpTwoScoreMatrix,
		TileSize:       32,
		StepSize:       32,
		GapPenalty:     -600,
		OpenGapPenalty: -10,
	}
	memoryPool := MatrixPoolMemory(defaultMatrixSize)
	m, trace := swMatrixSetup(10000)
	scoreTwoBit, alignmentPathTwoBit, targetEndTwoBit, queryEndTwoBit := TwoBitLocalRightAlign(dnaTwoBit.NewTwoBit(seqOneA), dnaTwoBit.NewTwoBit(seqOneB), config, memoryPool)
	scoreDnaBase, alignmentPathDnaBase, targetStartDnaBase, targetEndDnaBase, queryStartDnaBase, queryEndDnaBase := RightLocal(seqOneA, seqOneB, align.HumanChimpTwoScoreMatrix, config.GapPenalty, m, trace)

	expextedScore := int64(1244)
	expectedAlign := "1X15="
	expectedTargetStart := 0
	expectedTargetEnd := 16
	expectedQueryStart := 0
	expectedQueryEnd := 16

	t.Logf("score=%d, alignment=%s, refEnd=%d,  queryEnd=%d\n", scoreTwoBit, cigar.ByteCigarToString(alignmentPathTwoBit), targetEndTwoBit, queryEndTwoBit)

	if scoreTwoBit != expextedScore || scoreTwoBit != scoreDnaBase {
		t.Errorf("Error: Right direction expected score %d, got %d", expextedScore, scoreTwoBit)
	}
	if cigar.ByteCigarToString(alignmentPathTwoBit) != expectedAlign || cigar.ByteCigarToString(alignmentPathTwoBit) != cigar.ByteCigarToString(alignmentPathDnaBase) {
		t.Errorf("Error: Right direction expected score %s, got %s", cigar.ByteCigarToString(alignmentPathTwoBit), expectedAlign)
	}

	if expectedTargetStart != targetStartDnaBase {
		t.Errorf("Error: Right target start position %d, != %d", expectedTargetStart, targetStartDnaBase)
	}

	if targetEndTwoBit != expectedTargetEnd || targetEndTwoBit != targetEndDnaBase {
		t.Errorf("Error: Right target end position %d, != %d", expectedTargetEnd, targetEndTwoBit)
	}

	if expectedQueryStart != queryStartDnaBase {
		t.Errorf("Error: Right query start position %d, != %d", expectedQueryStart, queryStartDnaBase)
	}

	if queryEndTwoBit != expectedQueryEnd || queryEndTwoBit != queryEndDnaBase {
		t.Errorf("Error: Right query end position %d, != %d", expectedQueryEnd, queryEndTwoBit)
	}
}
