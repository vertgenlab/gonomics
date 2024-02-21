package genomeGraph

import (
	"reflect"
	"sync"
	"testing"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

func TestAlignGraphTraversal(t *testing.T) {
	// Setup a simple genome graph for testing
	genome := EmptyGraph()

	var n0, n1, n2 Node
	var e0, e1 Edge

	// Make Nodes
	n0 = Node{
		Id:  0,
		Seq: dna.StringToBases("TTG"),
	}

	n1 = Node{
		Id:  1,
		Seq: dna.StringToBases("AA"),
	}

	n2 = Node{
		Id:  2,
		Seq: dna.StringToBases("TTTTTT"),
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
	n1.Prev = append(n1.Prev, Edge{Dest: &n0, Prob: 1}) // Update to include reverse direction
	n2.Prev = append(n2.Prev, Edge{Dest: &n0, Prob: 1}) // Update to include reverse direction

	genome.Nodes = append(genome.Nodes, n0, n1, n2)

	// Create mock pool
	pool := &sync.Pool{
		New: func() interface{} {
			return &dnaPool{} // Assuming dnaPool is defined
		},
	}

	// Simulate a read
	read := dna.StringToBases("TGAAT")

	// Initialize matrix and scoreKeeper
	matrix := NewSwMatrix(defaultMatrixSize) // Example initialization

	sk := scoreKeeper{}                  // Initialize properly
	dynamicScore := dynamicScoreKeeper{} // Initialize properly

	// Perform traversals in both directions
	rightAlign, rightScore, rightTargetStart, rightQueryStart, rightPath := AlignGraphTraversal(&genome.Nodes[0], read, 0, []uint32{}, len(read), read, HumanChimpTwoScoreMatrix, &matrix, &sk, &dynamicScore, pool, rightTraversal)
	leftAlign, leftScore, leftTargetStart, leftQueryStart, leftPath := AlignGraphTraversal(&genome.Nodes[2], read, len(genome.Nodes[2].Seq), []uint32{}, len(read), read, HumanChimpTwoScoreMatrix, &matrix, &sk, &dynamicScore, pool, leftTraversal)
	// Expected results for RightDirection
	expectedRightScore := int64(460) // Simplified expected score
	expectedRightAlignment := []cigar.ByteCigar{{RunLen: 5, Op: 'M'}}
	expectedRightTargetStart := 5
	expectedRightQueryStart := 5

	// Expected results for LeftDirection
	expectedLeftScore := int64(460) // Simplified expected score, assuming a lesser match
	expectedLeftAlignment := []cigar.ByteCigar{{RunLen: 5, Op: 'M'}}
	expectedLeftTargetStart := 5 // Assuming the alignment starts one base into the sequence for simplicity
	expectedLeftQueryStart := 5
	expectedLeftPath, expectedRightPath := []uint32{2}, []uint32{}

	// Assert RightDirection results
	if rightScore != expectedRightScore {
		t.Errorf("RightDirection: expected score %d, got %d", expectedRightScore, rightScore)
	}

	if cigar.ByteCigarToString(rightAlign) != cigar.ByteCigarToString(expectedRightAlignment) {
		t.Errorf("RightDirection: expected alignment %+v, got %+v", expectedRightAlignment, rightAlign)
	}

	if rightTargetStart != expectedRightTargetStart || rightQueryStart != expectedRightQueryStart {
		t.Errorf("RightDirection: expected target and query start (%d, %d), got (%d, %d)", expectedRightTargetStart, expectedRightQueryStart, rightTargetStart, rightQueryStart)
	}

	if len(rightPath) != len(expectedRightPath) {
		t.Errorf("Expected path %+v, got %+v", expectedRightPath, rightPath)
	}

	if len(leftPath) != len(expectedLeftPath) {
		t.Errorf("Expected path %+v, got %+v", expectedLeftPath, leftPath)
	}
	// Assert LeftDirection results
	if leftScore != expectedLeftScore {
		t.Errorf("LeftDirection: expected score %d, got %d", expectedLeftScore, leftScore)
	}
	if !reflect.DeepEqual(leftAlign, expectedLeftAlignment) {
		t.Errorf("LeftDirection: expected alignment %+v, got %+v", expectedLeftAlignment, leftAlign)
	}

	if leftTargetStart != expectedLeftTargetStart || leftQueryStart != expectedLeftQueryStart {
		t.Errorf("LeftDirection: expected target and query start (%d, %d), got (%d, %d)", expectedLeftTargetStart, expectedLeftQueryStart, leftTargetStart, leftQueryStart)
	}
}

func TestDynamicAln(t *testing.T) {
	// Setup
	scores := [][]int64{
		{2, -1, -1, -1},
		{-1, 2, -1, -1},
		{-1, -1, 2, -1},
		{-1, -1, -1, 2},
	}
	gapPen := int64(-2)
	matrix := NewSwMatrix(defaultMatrixSize) // Assume a constructor for MatrixAln
	var dynamicScore dynamicScoreKeeper

	cases := []struct {
		name      string
		alpha     []dna.Base
		beta      []dna.Base
		direction byte
		expScore  int64
		expCigar  []cigar.ByteCigar
		expStartI int
		expStartJ int
	}{
		{
			name:      "Rightward Traversal",
			alpha:     dna.StringToBases("ACGT"),
			beta:      dna.StringToBases("ACGT"),
			direction: rightTraversal,
			expScore:  8,
			expCigar:  []cigar.ByteCigar{{RunLen: 4, Op: 'M'}},
			expStartI: 4,
			expStartJ: 4,
		},
		{
			name:      "Leftward Traversal",
			alpha:     dna.StringToBases("ACGT"),
			beta:      dna.StringToBases("GACG"),
			direction: leftTraversal,
			expScore:  4, // Assuming the alignment starts at the 'G' for simplicity
			expCigar:  []cigar.ByteCigar{{RunLen: 3, Op: 'M'}},
			expStartI: 3,
			expStartJ: 4,
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			score, cig, startI, startJ := DynamicAln(c.alpha, c.beta, scores, &matrix, gapPen, &dynamicScore, c.direction)
			if score != c.expScore {
				t.Errorf("Expected score %d, got %d", c.expScore, score)
			}
			if cigar.ByteCigarToString(cig) != cigar.ByteCigarToString(c.expCigar) {
				t.Errorf("Expected cigar %v, got %v", c.expCigar, cig)
			}
			if startI != c.expStartI || startJ != c.expStartJ {
				t.Errorf("Expected start indices (%d, %d), got (%d, %d)", c.expStartI, c.expStartJ, startI, startJ)
			}
		})
	}
}

func TestNodeSeqTraversal(t *testing.T) {
	// Create a test node
	node := Node{
		Id:  0,
		Seq: dna.StringToBases("ACGT"),
	}

	// Test leftTraversal
	extension := 2
	position := 2
	seq := dna.StringToBases("TT")
	direction := leftTraversal
	expectedResult := dna.StringToBases("ACTT")
	result := nodeSeqTraversal(&node, extension, position, seq, []dna.Base{}, direction)
	if !reflect.DeepEqual(result, expectedResult) {
		t.Errorf("LeftTraversal: expected %v, got %v", expectedResult, result)
	}
	// Test rightTraversal
	extension = 2
	position = 1
	seq = dna.StringToBases("TT")
	direction = rightTraversal
	expectedResult = dna.StringToBases("TTCG")
	result = nodeSeqTraversal(&node, extension, position, seq, []dna.Base{}, direction)
	if !reflect.DeepEqual(result, expectedResult) {
		t.Errorf("RightTraversal: expected %v, got %v", expectedResult, result)
	}
}
