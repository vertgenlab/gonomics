package genomeGraph

// func TestExtendToTheLeftHelper(t *testing.T) {
// 	// Create a test node
// 	type Edge struct {
// 		Node *Node
// 	}

// 	node := &Node{
// 		Id:        1,
// 		ColumnId:  0,
// 		Seq:       []dna.Base{},
// 		SeqTwoBit: &dna.TwoBitSequence{Len: 100},
// 		Prev:      []GenomeGraph.Edge{},
// 		Next:      []GenomeGraph.Edge{},

// 	}

// 	// Create a test read
// 	read := fastq.FastqBig{
// 		Rainbow:    []uint32{1, 2, 3, 4, 5},
// 		RainbowRc:  []uint32{5, 4, 3, 2, 1},
// 	}

// 	// Create a test nextPart
// 	nextPart := &SeedDev{
// 		QueryStart:   10,
// 		PosStrand:    true,
// 		TotalLength:  20,
// 	}

// 	// Call the function being tested
// 	result := extendToTheLeftHelper(node, read, nextPart)

// 	// Perform assertions on the result
// 	if len(result) != 1 {
// 		log.Fatalf("Expected 1 result, but got %d", len(result))
// 	}

// 	// Add more assertions as needed
// 	// ...

// 	// Add more test cases as needed
// 	// ...
// }
