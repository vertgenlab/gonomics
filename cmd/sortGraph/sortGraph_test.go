package main

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/genomeGraph"
	"io/ioutil"
	"os"
	"testing"
)

func TestGetSortOrder(t *testing.T) {
	sortedGraph := test(MakeTestGraph())
	correctOrder := []uint32{0, 1, 2, 3, 4}

	for i, node := range sortedGraph.Nodes {
		if node.Id != correctOrder[i] {
			t.Error("problem with sortGraph")
		}
	}
}

func test(graph *genomeGraph.GenomeGraph) (sortedGraph *genomeGraph.GenomeGraph) {
	input, err := ioutil.TempFile("", "sortGraphTest")
	exception.PanicOnErr(err)
	output, err := ioutil.TempFile("", "sortGraphTest")
	exception.PanicOnErr(err)

	input.Close()
	output.Close()

	genomeGraph.Write(input.Name(), graph)
	sortGraph(input.Name(), output.Name())
	answer := genomeGraph.Read(output.Name())

	err = os.Remove(input.Name())
	exception.PanicOnErr(err)
	err = os.Remove(output.Name())
	exception.PanicOnErr(err)
	return answer
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
func MakeDisContigTestGraph() *genomeGraph.GenomeGraph {
	graph := genomeGraph.EmptyGraph()

	var n0, n1, n2, n3, n4 genomeGraph.Node
	var f0, f2, f5, r0, r2, r5 genomeGraph.Edge

	// Make Nodes
	n0 = genomeGraph.Node{
		Id:  0,
		Seq: dna.StringToBases("ATG")}

	n1 = genomeGraph.Node{
		Id:  1,
		Seq: dna.StringToBases("CG")}

	n2 = genomeGraph.Node{
		Id:  2,
		Seq: dna.StringToBases("A")}

	n3 = genomeGraph.Node{
		Id:  3,
		Seq: dna.StringToBases("T")}

	n4 = genomeGraph.Node{
		Id:  4,
		Seq: dna.StringToBases("TAA")}

	// Make Edges
	f0 = genomeGraph.Edge{
		Dest: &n1,
		Prob: 1}

	f2 = genomeGraph.Edge{
		Dest: &n4,
		Prob: 1}

	f5 = genomeGraph.Edge{
		Dest: &n4,
		Prob: 1}

	r0 = genomeGraph.Edge{
		Dest: &n0,
		Prob: 1}

	r2 = genomeGraph.Edge{
		Dest: &n2,
		Prob: 1}

	r5 = genomeGraph.Edge{
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

// Sam Header:
//@HD     VN:1.6  SO:coordinate
//@SQ     SN:n0 LN:3
//@SQ     SN:n1 LN:2
//@SQ     SN:n2 LN:1
//@SQ     SN:n3 LN:1
//@SQ     SN:n4 LN:3

// Test Functions
func MakeTestGraph() *genomeGraph.GenomeGraph {
	graph := genomeGraph.EmptyGraph()

	var n0, n1, n2, n3, n4 genomeGraph.Node
	var e0, e1, e2, e3, e4, e5 genomeGraph.Edge

	// Make Nodes
	n0 = genomeGraph.Node{
		Id:  0,
		Seq: dna.StringToBases("ATG")}

	n1 = genomeGraph.Node{
		Id:  1,
		Seq: dna.StringToBases("CG")}

	n2 = genomeGraph.Node{
		Id:  2,
		Seq: dna.StringToBases("A")}

	n3 = genomeGraph.Node{
		Id:  3,
		Seq: dna.StringToBases("T")}

	n4 = genomeGraph.Node{
		Id:  4,
		Seq: dna.StringToBases("TAA")}

	// Make Edges
	e0 = genomeGraph.Edge{
		Dest: &n1,
		Prob: 1}

	e1 = genomeGraph.Edge{
		Dest: &n2,
		Prob: 0.05}

	e2 = genomeGraph.Edge{
		Dest: &n4,
		Prob: 1}

	e3 = genomeGraph.Edge{
		Dest: &n4,
		Prob: 0.8}

	e4 = genomeGraph.Edge{
		Dest: &n3,
		Prob: 0.15}

	e5 = genomeGraph.Edge{
		Dest: &n4,
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
