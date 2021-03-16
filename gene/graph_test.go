package gene

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/genomeGraph"
	"testing"
)

// TestGraph Structure
//             n2          e0 = 1
//         e1/    \e2      e1 = 0.05
//      e0  /  e3  \       e2 = 1
//  n0 --- n1 ----- n4     e3 = 0.8
//          \      /       e4 = 0.15
//         e4\    /e5      e5 = 1
//             n3
//
//               A
//             /    \
//            /      \
//  ATG --- CG ----- TAA
//            \      /
//             \    /
//               T

func MakeTestGraph() *genomeGraph.GenomeGraph {
	graph := genomeGraph.EmptyGraph()

	var n0, n1, n2, n3, n4 *genomeGraph.Node
	var e0, e1, e2, e3, e4, e5 genomeGraph.Edge

	// Make Nodes
	n0 = &genomeGraph.Node{
		Id:  0,
		Seq: dna.StringToBases("ATG")}

	n1 = &genomeGraph.Node{
		Id:  1,
		Seq: dna.StringToBases("CG")}

	n2 = &genomeGraph.Node{
		Id:  2,
		Seq: dna.StringToBases("A")}

	n3 = &genomeGraph.Node{
		Id:  3,
		Seq: dna.StringToBases("T")}

	n4 = &genomeGraph.Node{
		Id:  4,
		Seq: dna.StringToBases("TAA")}

	// Make Edges
	e0 = genomeGraph.Edge{
		Dest: n1,
		Prob: 1}

	e1 = genomeGraph.Edge{
		Dest: n2,
		Prob: 0.05}

	e2 = genomeGraph.Edge{
		Dest: n4,
		Prob: 1}

	e3 = genomeGraph.Edge{
		Dest: n4,
		Prob: 0.8}

	e4 = genomeGraph.Edge{
		Dest: n3,
		Prob: 0.15}

	e5 = genomeGraph.Edge{
		Dest: n4,
		Prob: 1}

	// Define Paths
	n0.Next = append(n0.Next, e0)
	n1.Next = append(n1.Next, e1, e3, e4)
	n1.Prev = append(n1.Prev, e0)
	n2.Next = append(n2.Next, e2)
	n2.Prev = append(n2.Prev, e1)
	n3.Next = append(n3.Next, e5)
	n3.Prev = append(n3.Prev, e4)
	n4.Prev = append(n4.Prev, e2, e3, e5)

	graph.Nodes = append(graph.Nodes, n0, n1, n2, n3, n4)

	for i := 0; i < len(graph.Nodes); i++ {
		graph.Nodes[i].SeqTwoBit = dnaTwoBit.NewTwoBit(graph.Nodes[i].Seq)
	}

	return graph
}

func TestSeqToPath(t *testing.T) {
	g := MakeTestGraph()

	// test 1
	query := dna.StringToBases("ATGCGTAA")
	ids, path, err := SeqToPath(query, g.Nodes[0], 0, *g)
	if err != nil {
		t.Error(err)
	}
	testSeq := extractSeq(path.Nodes[ids[0]])
	expectedRev := dna.StringToBases("TAACGATG")
	actualRev := getRevNodeSeq(path.Nodes[ids[0]])

	if dna.CompareSeqsIgnoreCase(query, testSeq) != 0 ||
		dna.CompareSeqsIgnoreCase(expectedRev, actualRev) != 0 {
		t.Error("problem with SeqToPath")
	}

	for _, id := range ids {
		if path.Nodes[id].Id != g.Nodes[id].Id {
			t.Error("problem with SeqToPath")
		}
	}

	if countNonNil(path) != len(ids) {
		t.Error("problem with SeqToPath")
	}

	// test 2
	query = dna.StringToBases("TGCG")
	ids, path, err = SeqToPath(query, g.Nodes[0], 1, *g)
	if err != nil {
		t.Error(err)
	}
	testSeq = extractSeq(path.Nodes[ids[0]])

	expectedRev = dna.StringToBases("CGTG")
	actualRev = getRevNodeSeq(path.Nodes[ids[0]])

	if dna.CompareSeqsIgnoreCase(query, testSeq) != 0 ||
		dna.CompareSeqsIgnoreCase(expectedRev, actualRev) != 0 {
		t.Error("problem with SeqToPath")
	}

	for _, id := range ids {
		if path.Nodes[id].Id != g.Nodes[id].Id {
			t.Error("problem with SeqToPath")
		}
	}

	if countNonNil(path) != len(ids) {
		t.Error("problem with SeqToPath")
	}

	// test 3
	query = dna.StringToBases("ATGCGT")
	ids, path, err = SeqToPath(query, g.Nodes[0], 0, *g)
	if err != nil {
		t.Error(err)
	}
	testSeq = extractSeq(path.Nodes[ids[0]])

	expectedRev = dna.StringToBases("TCGATG")
	actualRev = getRevNodeSeq(path.Nodes[ids[0]])

	if dna.CompareSeqsIgnoreCase(query, testSeq) != 0 ||
		dna.CompareSeqsIgnoreCase(expectedRev, actualRev) != 0 {
		t.Error("problem with SeqToPath")
	}

	for _, id := range ids {
		if path.Nodes[id].Id != g.Nodes[id].Id {
			t.Error("problem with SeqToPath")
		}
	}

	if countNonNil(path) != len(ids) {
		t.Error("problem with SeqToPath")
	}
}

func extractSeq(start *genomeGraph.Node) []dna.Base {
	var answer []dna.Base
	var currNode *genomeGraph.Node = start
	answer = append(answer, start.Seq...)

	for currNode.Next != nil {
		answer = append(answer, currNode.Next[0].Dest.Seq...)
		currNode = currNode.Next[0].Dest
	}

	return answer
}

func getRevNodeSeq(start *genomeGraph.Node) []dna.Base {
	currNode := start
	for currNode.Next != nil { // move to end
		currNode = currNode.Next[0].Dest
	}
	var answer []dna.Base
	for currNode.Prev != nil {
		answer = append(answer, currNode.Seq...)
		currNode = currNode.Prev[0].Dest
	}
	answer = append(answer, currNode.Seq...)
	return answer
}

func countNonNil(g genomeGraph.GenomeGraph) int {
	var answer int
	for i := range g.Nodes {
		if g.Nodes[i] != nil {
			answer++
		}
	}
	return answer
}
