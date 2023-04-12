package genomeGraph

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
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

	var discontiguousGraph *GenomeGraph = MakeDisContigTestGraph()
	brokenGraph := BreakNonContiguousGraph(discontiguousGraph.Nodes)
	if len(brokenGraph) != 2 {
		t.Errorf("ERROR: Problem with breaking up discontigous graph into subgraphs")
	}

}

// TestGraph Structure
//             n2          e0 = 1
//         e1     \e2
//      e0     e3  \       e2 = 1
//  n0 --- n1       n4
//                 /
//         e4     /e5      e5 = 1
//             n3
//
//               A
//                  \
//                   \
//  ATG --- CG 		 TAA
//                   /
//                  /
//               T
func MakeDisContigTestGraph() *GenomeGraph {
	graph := EmptyGraph()

	var n0, n1, n2, n3, n4 Node
	var f0, f2, f5, r0, r2, r5 Edge

	// Make Nodes
	n0 = Node{
		Id:  0,
		Seq: dna.StringToBases("ATG")}

	n1 = Node{
		Id:  1,
		Seq: dna.StringToBases("CG")}

	n2 = Node{
		Id:  2,
		Seq: dna.StringToBases("A")}

	n3 = Node{
		Id:  3,
		Seq: dna.StringToBases("T")}

	n4 = Node{
		Id:  4,
		Seq: dna.StringToBases("TAA")}

	// Make Edges
	f0 = Edge{
		Dest: &n1,
		Prob: 1}

	f2 = Edge{
		Dest: &n4,
		Prob: 1}

	f5 = Edge{
		Dest: &n4,
		Prob: 1}

	r0 = Edge{
		Dest: &n0,
		Prob: 1}

	r2 = Edge{
		Dest: &n2,
		Prob: 1}

	r5 = Edge{
		Dest: &n3,
		Prob: 1}

	// Define Paths
	n0.Next = append(n0.Next, f0)
	n1.Prev = append(n1.Prev, r0)
	n2.Next = append(n2.Next, f2)
	n3.Next = append(n3.Next, f5)
	n4.Prev = append(n4.Prev, r2, r5)

	graph.Nodes = append(graph.Nodes, n0, n1, n2, n3, n4)

	for i := 0; i < len(graph.Nodes); i++ {
		graph.Nodes[i].SeqTwoBit = dnaTwoBit.NewTwoBit(graph.Nodes[i].Seq)
	}

	return graph
}
