package simpleGraph

import (
	"testing"
)

func TestGetSortOrder(t *testing.T) {
	graph := MakeTestGraph()

	// Just getting two sets of the same nodes to make sure they will separate
	nodesToSort := append(graph.Nodes)
	order := breadthFirstSearch(nodesToSort)

	correctOrder := []uint32{0, 1, 2, 3, 4}

	for i := 0; i < len(order); i++ {
		if order[i] == correctOrder[i] {
			continue
		}
		t.Errorf("Error: Problem with graph sort \n Expected Order : %v \n Received Order : %v", correctOrder, order)
	}
}
