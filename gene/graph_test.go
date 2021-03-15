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

// TODO: check prevs
func TestSeqToPath(t *testing.T) {
	g := MakeTestGraph()

	// test 1
	start, _, err := SeqToPath(dna.StringToBases("ATGCGTAA"), g.Nodes[0], 0, *g)
	if err != nil {
		t.Error(err)
	}

	expected := dna.StringToBases("ATGCGTAA")
	testSeq := extractSeq(start)

	if dna.CompareSeqsIgnoreCase(expected, testSeq) != 0 {
		t.Error("problem with SeqToPath")
	}

	// test 2
	start, _, err = SeqToPath(dna.StringToBases("TGCG"), g.Nodes[0], 1, *g)
	if err != nil {
		t.Error(err)
	}

	expected = dna.StringToBases("TGCG")
	testSeq = extractSeq(start)

	if dna.CompareSeqsIgnoreCase(expected, testSeq) != 0 {
		t.Error("problem with SeqToPath")
	}

	// test 3
	start, _, err = SeqToPath(dna.StringToBases("ATGCGT"), g.Nodes[0], 0, *g)
	if err != nil {
		t.Error(err)
	}

	expected = dna.StringToBases("ATGCGT")
	testSeq = extractSeq(start)

	if dna.CompareSeqsIgnoreCase(expected, testSeq) != 0 {
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
