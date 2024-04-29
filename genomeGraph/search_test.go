package genomeGraph

import (
	"testing"
	"time"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
)

// TestGetRightBases tests the getRightBases function for correct behavior.
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
	rightTargetBases := getRightBases(n1, 10, 2, rightSeq, dna.StringToBases(""))

	if dna.BasesToString(rightTargetBases) != dna.BasesToString(rightExpected) {
		t.Errorf("getRightBases() = %v, want %v", dna.BasesToString(rightTargetBases), dna.BasesToString(rightExpected))
	}
	leftSeq := dna.StringToBases("A")
	leftExpected := dna.StringToBases("CGATCGA")
	leftTargetBases := getLeftBases(n2, 7, 7, leftSeq, dna.StringToBases(""))

	if dna.BasesToString(leftExpected) != dna.BasesToString(leftTargetBases) {
		t.Errorf("getLeftTargetBases() = %v, want %v", dna.BasesToString(leftTargetBases), dna.BasesToString(leftExpected))
	}
}

func TestRightDynamicAln(t *testing.T) {
	alpha := dna.StringToBases("ACGTA")
	beta := dna.StringToBases("ACGATCGAT")

	settings := &GraphSettings{
		ScoreMatrix: align.DefaultScoreMatrix,
		GapPenalty:  -100,
		TileSize:    32,
		StepSize:    32,
		Extension:   150,
	}

	expectedScore := int64(291)
	expectedCigar := []cigar.ByteCigar{{RunLen: 3, Op: 'M'}}
	expectedEndAlpha := 3
	expectedEndBeta := 3

	memory := MatrixPoolMemory(defaultMatrixSize)
	cache := memory.Get().(*MatrixMemoryPool)
	defer memory.Put(cache)

	score, cig, endAlpha, endBeta := RightDynamicAln(alpha, beta, settings, memory)

	if score != expectedScore {
		t.Errorf("Expected score: %d, got: %d", expectedScore, score)
	}

	if cigar.ByteCigarToString(cig) != cigar.ByteCigarToString(expectedCigar) {
		t.Errorf("Error: expected CIGAR: %v, got: %v", expectedCigar, cig)
	}

	if endAlpha != expectedEndAlpha {
		t.Errorf("Error: expected endAlpha: %d, got: %d", expectedEndAlpha, endAlpha)
	}

	if endBeta != expectedEndBeta {
		t.Errorf("Error: expected endBeta: %d, got: %d", expectedEndBeta, endBeta)
	}
}

func TestAlignGraphTraversal(t *testing.T) {
	// Setup a simple genome graph for testing
	genome := EmptyGraph()

	var n0, n1, n2 Node
	var e0, e1 Edge
	/*
		Graph nodes: 		AAAAAT -> ATCG -> CGTTTTTT
		Input Read: 		ATATCGCGT
		Output Expected:	AT->ATCG->CGT
	*/
	n0 = Node{
		Id:  0,
		Seq: dna.StringToBases("AAAAAT"),
	}

	n1 = Node{
		Id:  1,
		Seq: dna.StringToBases("AAAATCGATAGGGGAAAATTT"),
	}

	n2 = Node{
		Id:  2,
		Seq: dna.StringToBases("CGTTTTTT"),
	}

	// Make Edges
	e0 = Edge{
		Dest: &n1,
		Prob: 1,
	}

	e1 = Edge{
		Dest: &n2,
		Prob: 1,
	}
	// Define Paths
	n0.Next = append(n0.Next, e0, e1)
	n1.Prev = append(n1.Prev, Edge{Dest: &n0, Prob: 1})
	n2.Prev = append(n2.Prev, Edge{Dest: &n0, Prob: 1})

	genome.Nodes = append(genome.Nodes, n0, n1, n2)

	memory := MatrixPoolMemory(defaultMatrixSize)
	cache := memory.Get().(*MatrixMemoryPool)
	defer memory.Put(cache)

	read := dna.StringToBases("ATATCGCGT")

	scoreKeeper := scoreKeeper{}

	// expectedRightPath :=  []uint32{1,2}
	// Expected results for RightDirection
	expectedRightScore := int64(850) // Simplified expected score
	expectedRightAlignment := []cigar.ByteCigar{{RunLen: 9, Op: 'M'}}
	expectedRightTargetEnd := 9
	expectedRightQueryEnd := 9

	expectedLeftScore := int64(850) // Simplified expected score
	expectedLeftAlignment := []cigar.ByteCigar{{RunLen: 9, Op: 'M'}}
	expectedLeftTargetStart := 2
	expectedLeftQueryStart := 0

	settings := &GraphSettings{
		ScoreMatrix: align.HumanChimpTwoScoreMatrix,
		GapPenalty:  -100,
		TileSize:    32,
		Extension:   len(read),
	}

	rightAlign, rightScore, rightTargetEnd, rightQueryEnd, _ := RightAlignTraversal(&genome.Nodes[0], read, 0, []uint32{}, read, settings, &scoreKeeper, memory)
	leftAlign, leftScore, leftTargetStart, leftQueryStart, _ := LeftAlignTraversal(&genome.Nodes[0], read, 20, []uint32{}, read, settings, scoreKeeper, memory)
	if rightScore != expectedRightScore {
		t.Errorf("Error: Right direction expected score %d, got %d", expectedRightScore, rightScore)
	}
	if leftScore != expectedLeftScore {
		t.Errorf("Error: Left direction expected score %d, got %d", expectedLeftScore, leftScore)
	}
	if cigar.ByteCigarToString(rightAlign) != cigar.ByteCigarToString(expectedRightAlignment) {
		t.Errorf("Error: Right direction expected alignment %+v, got %+v", expectedRightAlignment, rightAlign)
	}
	if cigar.ByteCigarToString(leftAlign) != cigar.ByteCigarToString(expectedLeftAlignment) {
		t.Errorf("Error: Left direction expected alignment %+v, got %+v", expectedLeftAlignment, leftAlign)
	}
	if expectedRightTargetEnd != rightTargetEnd || rightQueryEnd != expectedRightQueryEnd {
		t.Errorf("Error: Right direction: expected target and query start (%d, %d), got (%d, %d)", expectedRightTargetEnd, rightQueryEnd, rightTargetEnd, rightQueryEnd)
	}
	if expectedLeftTargetStart != leftTargetStart || leftQueryStart != expectedLeftQueryStart {
		t.Errorf("Error: Left direction expected target and query start (%d, %d), got (%d, %d)", expectedLeftTargetStart, leftQueryStart, leftTargetStart, leftQueryStart)
	}

}

func BenchmarkGirafAlignment(b *testing.B) {

	var tileSize int = 32
	var stepSize int = 32
	config := &GraphSettings{
		ScoreMatrix:    align.HumanChimpTwoScoreMatrix,
		GapPenalty:     -600,
		OpenGapPenalty: -150,
		TileSize:       tileSize,
		StepSize:       stepSize,
	}
	config.MaxMatch, config.MinMatch, config.LeastSevereMismatch, config.LeastSevereMatchMismatchChange = MismatchStats(config.ScoreMatrix)
	genome := Read("testdata/bigGenome.sg")
	fqOne := "testdata/simReads_R1.fq"
	tiles := IndexGenomeIntoMap(genome.Nodes, tileSize, stepSize)

	fqs := []fastq.FastqBig{}
	file := fileio.NewByteReader(fqOne)
	for read, done := fastq.ReadFqBig(file); !done; read, done = fastq.ReadFqBig(file) {
		fqs = append(fqs, read)
	}

	seedPool := NewMemSeedPool()
	seedBuildHelper := newSeedBuilder()
	scorekeeper := scoreKeeper{}

	memory := MatrixPoolMemory(defaultMatrixSize)
	cache := memory.Get().(*MatrixMemoryPool)
	defer memory.Put(cache)

	results := []*giraf.Giraf{}
	var correct int = 0
	var start, stop time.Time
	for n := 0; n < b.N; n++ {

		start = time.Now()
		for read := 0; read < len(fqs); read++ {
			giraf := GraphSmithWatermanToGiraf(genome, fqs[read], tiles, config, memory, &seedPool, scorekeeper, seedBuildHelper)
			if IsCorrectCoord(giraf) {
				correct++
			}
			results = append(results, giraf)
		}
		stop = time.Now()
		duration := stop.Sub(start)
		b.Logf("Aligned %v out of %v correctly...\n", correct, len(results))
		b.Logf("Aligned %d reads in %s (%.1f reads per second).\n", len(results), duration, float64(len(fqs))/duration.Seconds())
		b.ReportAllocs()
	}
}
