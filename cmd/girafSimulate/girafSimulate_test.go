package main

import (
	"io/ioutil"
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/genomeGraph"
	"github.com/vertgenlab/gonomics/giraf"
)

func TestGirafSimulate(t *testing.T) {
	graph := MakeTestGraph()
	out, err := ioutil.TempFile("", "GirafSimulate_test")
	exception.PanicOnErr(err)
	out.Close()
	girafSimulate(graph, 10, 3, 0, 1, 0.2, out.Name(), false)

	reads := giraf.Read(out.Name())
	if len(reads) != 10 {
		t.Error("problem with girafSimulate")
	}
	for i := range reads {
		if len(reads[i].Seq) != 3 {
			t.Error("problem with girafSimulate")
		}
	}
	fileio.EasyRemove(out.Name())
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
