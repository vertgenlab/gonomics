package genomeGraph

import (
	"reflect"
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers/parse"
)

func TestLeftDynamicAln(t *testing.T) {
	alpha := dna.StringToBases("ATGC")
	beta := dna.StringToBases("ATCC")
	scores := align.HumanChimpTwoScoreMatrix
	gapPen := int64(-10)

	expectedScore := int64(260)
	expectedCigar := cigar.ByteCigarToString([]cigar.ByteCigar{{RunLen: 1, Op: cigar.Match}, {RunLen: 1, Op: cigar.Insertion}, {RunLen: 1, Op: cigar.Deletion}, {RunLen: 2, Op: cigar.Match}})
	expectedI := 0
	expectedJ := 0

	// Setup
	matrix := NewSwMatrix((len(alpha) + len(beta)))
	dynamicScore := dynamicScoreKeeper{}

	// Call the function under test
	score, cig, i, j := LeftDynamicAln(alpha, beta, scores, &matrix, gapPen, dynamicScore)

	// Assertions
	if score != expectedScore {
		t.Errorf("Error: Expected score %d, but got %d", expectedScore, score)
	}
	result := cigar.ByteCigarToString(cig)
	if expectedCigar != result {
		t.Errorf("Error: Expected cigar %v, but got %v", expectedCigar, result)
	}
	if i != expectedI || j != expectedJ {
		t.Errorf("Error: Expected i, j = %d, %d, but got %d, %d", expectedI, expectedJ, i, j)
	}
}

func TestRightAlignTraversal(t *testing.T) {
	// Test Data
	graph := EmptyGraph()
	query := dna.StringToBases("ACTGACTGGGTTTTGGGAT") // Simple query sequence

	var n0, n1, n2 Node

	// Make Nodes
	n0 = Node{
		Id:  0,
		Seq: dna.StringToBases("ATGACTGACTGG")}

	n1 = Node{
		Id:  1,
		Seq: dna.StringToBases("GTTTTGG")}

	n2 = Node{
		Id:  2,
		Seq: dna.StringToBases("ATA")}

	AddEdge(&n0, &n1, 1)
	AddEdge(&n1, &n2, 1)

	graph.Nodes = []Node{n0, n1, n2}
	config := &GraphSettings{
		Extention: 40, // Adjust for test if needed
		ScoreMatrix: [][]int64{
			{5, -1, -1, -1},
			{-1, 5, -1, -1},
			{-1, -1, 5, -1},
			{-1, -1, -1, 5},
		},
		GapPenalty: -1,
	}

	// Expected Results
	expectedScore := int64(35)
	expectedCigar := cigar.ByteCigarToString([]cigar.ByteCigar{
		{RunLen: 9, Op: cigar.Insertion},
		{RunLen: 5, Op: cigar.Match},
		{RunLen: 1, Op: cigar.Insertion},
		{RunLen: 4, Op: cigar.Match},
	})
	expectedTargetEnd := 9
	expectedQueryEnd := 19
	expectedPath1 := []uint32{} // Path for node1 traversal
	expectedPath2 := []uint32{}

	// Test with node1 (traversal)
	matrix := NewSwMatrix(defaultMatrixSize)
	sk := scoreKeeper{}
	dynamicScore := dynamicScoreKeeper{}
	dnaPool := NewDnaPool()

	cigar1, score1, targetEnd1, queryEnd1, path1 := RightAlignTraversal(&n1, []dna.Base{}, 0, []uint32{}, query, config, &matrix, sk, dynamicScore, &dnaPool)

	if score1 != expectedScore {
		t.Errorf("Error: Expected score %d for node1, got %d", expectedScore, score1)
	}
	result1 := cigar.ByteCigarToString(cigar1)
	if cigar.ByteCigarToString(cigar1) != expectedCigar {
		t.Errorf("Error: Expected cigar %s, but got %s", expectedCigar, result1)
	}
	if targetEnd1 != expectedTargetEnd || queryEnd1 != expectedQueryEnd {
		t.Errorf("Error: Expected targetEnd, queryEnd = %d, %d for node1, got %d, %d", expectedTargetEnd, expectedQueryEnd, targetEnd1, queryEnd1)
	}
	if !reflect.DeepEqual(path1, expectedPath1) {
		t.Errorf("Error: Expected path %v for node1, got %v", expectedPath1, path1)
	}

	// Test with node2 (no traversal, direct alignment)
	cigar2, score2, targetEnd2, queryEnd2, path2 := RightAlignTraversal(&n1, []dna.Base{}, 0, []uint32{}, query, config, &matrix, sk, dynamicScore, &dnaPool)

	// Similar assertions for node2
	if score2 != expectedScore {
		t.Errorf("Error: Expected score %d for node2, got %d", expectedScore, score2)
	}
	result2 := cigar.ByteCigarToString(cigar2)
	if result2 != expectedCigar {
		t.Errorf("Error: Expected cigar %s, but got %s", expectedCigar, result1)
	}

	if targetEnd2 != expectedTargetEnd || queryEnd2 != expectedQueryEnd {
		t.Errorf("Error: Expected targetEnd, queryEnd = %d, %d, but got %d, %d", targetEnd2, queryEnd2, expectedTargetEnd, expectedQueryEnd)
	}
	if !reflect.DeepEqual(path2, expectedPath2) {
		t.Errorf("Error: Expected path RightTraversal %v for node1, got %v", expectedPath2, path2)
	}
}

func TestRightDynamicAln(t *testing.T) {
	// Test Case
	alpha := dna.StringToBases("GATTACA")
	beta := dna.StringToBases("GATACCA")
	config := &GraphSettings{
		ScoreMatrix: [][]int64{
			{2, -1, -1, -1},
			{-1, 2, -1, -1},
			{-1, -1, 2, -1},
			{-1, -1, -1, 2},
		},
		GapPenalty: -3, // Adjust as needed
	}

	expectedScore := int64(8) // Based on the scoring matrix and gap penalty
	expectedCigar := cigar.ByteCigarToString([]cigar.ByteCigar{{RunLen: 7, Op: cigar.Match}})
	expectedI := 7
	expectedJ := 7

	matrix := NewSwMatrix(len(alpha) + len(beta))
	dynamicScore := dynamicScoreKeeper{}

	// Call the function under test
	score, cig, i, j := RightDynamicAln(alpha, beta, config, &matrix, dynamicScore)

	// Assertions
	if score != expectedScore {
		t.Errorf("Error: Expected score %d, but got %d", expectedScore, score)
	}
	result := cigar.ByteCigarToString(cig)
	if expectedCigar != result {
		t.Errorf("Error: Expected cigar %v, but got %v", expectedCigar, result)
	}
	if i != expectedI || j != expectedJ {
		t.Errorf("Error: Expected i, j = %d, %d, but got %d, %d", expectedI, expectedJ, i, j)
	}
}

func BenchmarkGirafAlignment(b *testing.B) {
	// Setup (outside of benchmark timing)
	var numberOfReads int = 10
	var readLength int = 150
	var mutations int = 1

	config := &GraphSettings{
		ScoreMatrix:    align.HumanChimpTwoScoreMatrix,
		GapPenalty:     -600,
		OpenGapPenalty: -150,
		TileSize:       32,
		StepSize:       32,
	}

	genome := Read("testdata/bigGenome.sg")
	tiles := IndexGenomeIntoMap(genome.Nodes, config.TileSize, config.StepSize)

	simReads := RandomPairedReads(genome, readLength, numberOfReads, mutations)
	fqOne := "testdata/simReads_R1.fq"
	fqTwo := "testdata/simReads_R2.fq"

	fastq.WritePair(fqOne, fqTwo, simReads)

	// Pre-load reads to avoid file I/O inside the benchmark
	fqs := []fastq.FastqBig{}
	file := fileio.NewByteReader(fqOne)
	for read, done := fastq.ReadFqBig(file); !done; read, done = fastq.ReadFqBig(file) {
		fqs = append(fqs, read)
	}
	defer file.Close()

	// Reusable objects (created once, used throughout)
	matrix := NewSwMatrix(defaultMatrixSize)
	seedPool := NewMemSeedPool()
	dnaPool := NewDnaPool()
	seedBuildHelper := newSeedBuilder()

	b.ResetTimer()

	var correct int
	for n := 0; n < b.N; n++ {
		correct = 0
		for i := 0; i < len(fqs); i++ {
			scorekeeper := scoreKeeper{} // Create new scorekeepers per alignment
			dynamicKeeper := dynamicScoreKeeper{}
			result := GraphSmithWatermanToGiraf(genome, fqs[i], tiles, config, &matrix, &seedPool, &dnaPool, scorekeeper, dynamicKeeper, seedBuildHelper)
			if isCorrectCoord(result) {
				correct++
			}
		}
	}
	// Stop benchmark timer & report stats
	b.StopTimer()
	// Report additional stats
	b.ReportMetric(float64(correct)/float64(len(fqs)), "Accuracy")
	b.ReportMetric(float64(len(simReads))/b.Elapsed().Seconds(), "ReadsPerSecond")
	b.Logf("Mapping results: %d / %d reads correctly", correct, len(fqs))

	// Cleanup
	fileio.EasyRemove("testdata/simReads_R1.fq")
	fileio.EasyRemove("testdata/simReads_R2.fq")
}

func isCorrectCoord(result *giraf.Giraf) bool {
	name := strings.Split(result.QName, "_")
	// TODO: Need better logic to look at node
	return parse.StringToInt(name[1]) == result.Path.TStart-result.QStart && parse.StringToInt(name[3]) == result.Path.TEnd && cigar.QueryRunLen(result.Cigar) == 150
}
