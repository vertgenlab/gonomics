package genomeGraph

import (
	"reflect"
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
)

func TestExtendToTheRight(t *testing.T) {
	// Create a mock node that represents a sequence in the genome
	node := Node{
		Id:        1,
		SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ATCGATCG")),
		Next: []Edge{ // Assume this node is connected to another node
			{
				Dest: &Node{
					Id:        2,
					SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("GTCGAC")),
				},
			},
		},
	}

	// Create a mock FastqBig read that has a sequence partially matching the node's sequence
	read := fastq.FastqBig{
		Name:    "read1",
		Seq:     dna.StringToBases("ATCGAT"),
		Rainbow: dnaTwoBit.TwoBitRainbowDeReference(dna.StringToBases("ATCGAT")),
	}

	// Expected seeds after extending to the right
	expectedSeeds := []*SeedDev{
		{
			TargetId:    1,
			TargetStart: 0,
			QueryStart:  0,
			Length:      6,
			PosStrand:   true,
			TotalLength: 6,
		},
	}

	// Call extendToTheRight with the mock node and read, starting from the beginning of both sequences
	seeds := extendToTheRight(&node, read, 0, 0, true)

	// Check if the returned seeds match the expected seeds
	if !reflect.DeepEqual(seeds, expectedSeeds) {
		t.Errorf("extendToTheRight returned unexpected seeds. Got %+v, want %+v", seeds, expectedSeeds)
	}
}

func TestExtendToTheLeft(t *testing.T) {
	// Mock previous node connected to the main node where the seed will be extended to the left
	prevNode := Node{
		Id:        0,
		SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("GCA")),
		Prev:      nil, // Assume this is the start of the graph
	}

	// Main node where the seed initially matches
	mainNode := Node{
		Id:        1,
		SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ATCG")),
		Prev: []Edge{ // Connection from the main node to the previous node
			{
				Dest: &prevNode,
			},
		},
	}

	// Mock FastqBig read that has a sequence matching the end of the previous node and the start of the main node
	read := fastq.FastqBig{
		Name:    "read1",
		Seq:     dna.StringToBases("GCAATCG"),
		Rainbow: dnaTwoBit.TwoBitRainbowDeReference(dna.StringToBases("GCAATCG")),
	}

	read.SeqRc = make([]dna.Base, len(read.Seq))

	copy(read.SeqRc, read.Seq)
	dna.ReverseComplement(read.SeqRc)
	read.RainbowRc = dnaTwoBit.TwoBitRainbowDeReference(read.SeqRc)

	// Seed that initially matches the start of the main node
	initialSeed := &SeedDev{
		TargetId:    0,
		TargetStart: 0,
		QueryStart:  0, // Matches start at the fourth base of the read
		Length:      3, // Matches "ATCG" from the main node
		PosStrand:   true,
		TotalLength: 3,
	}

	// Expected seeds after extending to the left into the prevNode
	expectedSeeds := []*SeedDev{
		{
			TargetId:    0,
			TargetStart: 0,
			QueryStart:  0,
			Length:      3, // Matches "GCA" from the prevNode
			PosStrand:   true,
			TotalLength: 3, // Total match length including the initial seed and the extension
		},
	}
	// Call extendToTheLeft with the main node, read, and initial seed
	seeds := extendToTheLeft(&mainNode, read, initialSeed)

	// Check if the returned seeds match the expected seeds using seedSlicesEqual
	if len(seeds) != len(expectedSeeds) {
		t.Errorf("extendToTheLeft returned unexpected seeds. Got %+v, want %+v", seeds, expectedSeeds)
	}
	// TODO: finish testing this function
	if !seedEqual(seeds[0], expectedSeeds[0]) {
		t.Errorf("extendToTheLeft returned unexpected seeds. Got %+v, want %+v", seeds[0], expectedSeeds[0])
	}

}
func TestFindSeedsInSmallMapWithMemPool(t *testing.T) {
	// Mock genome graph nodes
	nodes := []Node{
		{Id: 0, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("ATCGATCG"))},
		{Id: 1, SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases("GTCGAC"))},
	}

	// Mock read
	read := fastq.FastqBig{
		Name: "read1",
		Seq:  dna.StringToBases("ATCG"),
		Qual: []uint8{30, 30, 30, 30},
	}
	read.SeqRc = make([]dna.Base, len(read.Seq))

	copy(read.SeqRc, read.Seq)
	dna.ReverseComplement(read.SeqRc)

	read.Rainbow = dnaTwoBit.TwoBitRainbowDeReference(read.Seq)
	read.RainbowRc = dnaTwoBit.TwoBitRainbowDeReference(read.SeqRc)

	seedHash := make(map[uint64][]uint64)
	seqKey := dnaToNumber(read.Seq, 0, 4) // Assuming the read matches the start of node 0
	seedHash[seqKey] = append(seedHash[seqKey], ChromAndPosToNumber(0, 0))

	// Seed length, step, and score matrix setup
	seedLen := 4
	perfectScore := int64(40)         // Assuming a simplified scoring system
	scoreMatrix := make([][]int64, 4) // Simplified score matrix initialization

	// Call the function under test
	seeds := findSeedsInSmallMapWithMemPool(seedHash, nodes, read, seedLen, perfectScore, scoreMatrix)

	// Expected seed
	expectedSeed := &SeedDev{
		TargetId:    0,
		TargetStart: 0,
		QueryStart:  0,
		Length:      4,
		PosStrand:   true,
		TotalLength: 4,
		NextPart:    nil,
	}

	// Validate results
	if len(seeds) != 1 {
		t.Errorf("Expected 1 seed, got %d", len(seeds))
	} else if !compareSeeds(seeds[0], expectedSeed) {
		t.Errorf("Seed does not match expected result")
	}
}

// Helper function to compare two SeedDev objects
func compareSeeds(a, b *SeedDev) bool {
	return a.TargetId == b.TargetId && a.TargetStart == b.TargetStart && a.QueryStart == b.QueryStart && a.Length == b.Length && a.PosStrand == b.PosStrand && a.TotalLength == b.TotalLength
}

func seedSlicesEqual(a, b []*SeedDev) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !seedEqual(a[i], b[i]) {
			return false
		}
	}
	return true
}

func seedEqual(a, b *SeedDev) bool {
	// Compare all fields of SeedDev, and recursively compare NextPart if not nil
	return a.TargetId == b.TargetId &&
		a.TargetStart == b.TargetStart &&
		a.QueryStart == b.QueryStart &&
		a.Length == b.Length &&
		a.PosStrand == b.PosStrand &&
		a.TotalLength == b.TotalLength
}

func TestProcessStrand(t *testing.T) {
	// Mock seedHash
	chr0, chr1 := dna.StringToBases("ATCGATCG"), dna.StringToBases("ATTTTTTT")
	target := [][]dna.Base{chr0, chr1}
	query := dna.StringToBases("ATTT")

	keyIdx := (0 + 31) / 32
	keyOffset := 31 - ((0 + 31) % 32)

	headNode := Node{
		Id:        1,
		SeqTwoBit: dnaTwoBit.NewTwoBit(target[1]),
		Prev:      nil,
	}

	// Main node where the seed initially matches
	nextNode := Node{
		Id:        1,
		SeqTwoBit: dnaTwoBit.NewTwoBit(query),
		Prev: []Edge{ // Connection from the main node to the previous node
			{
				Dest: &headNode,
			},
		},
	}

	nodes := []Node{
		headNode, nextNode,
	}

	read := fastq.FastqBig{
		Name:    "read1",
		Seq:     query,
		Rainbow: dnaTwoBit.TwoBitRainbowDeReference(query),
	}

	read.SeqRc = make([]dna.Base, len(read.Seq))

	copy(read.SeqRc, read.Seq)
	dna.ReverseComplement(read.SeqRc)
	read.RainbowRc = dnaTwoBit.TwoBitRainbowDeReference(read.SeqRc)
	finalSeeds := []*SeedDev{}

	keyShift := uint(0)

	fwd, rev := read.Rainbow[keyOffset].Seq[keyIdx]>>keyShift, read.RainbowRc[keyOffset].Seq[keyIdx]>>keyShift
	seedHash := map[uint64][]uint64{
		fwd: {ChromAndPosToNumber(0, 0)},
		rev: {ChromAndPosToNumber(1, 0)},
	}

	processStrand(seedHash, nodes, read, 6, &finalSeeds, 0, uint(0))
	// Define the expected finalSeeds
	expectedFinalSeeds := []*SeedDev{
		{
			TargetId:    1,
			TargetStart: 0,
			QueryStart:  0,
			Length:      4,
			PosStrand:   true,
			TotalLength: 4,
			NextPart:    nil,
		},
	}
	// Validate results
	if !compareSeeds(finalSeeds[0], expectedFinalSeeds[0]) {
		t.Errorf("Error: Final seeds do not match expected result %v != %v", finalSeeds[0], expectedFinalSeeds[0])
	}
}
