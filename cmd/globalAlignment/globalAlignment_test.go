package main

import (
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"testing"
)

//initiate fake fasta files for testing the graph function
var toad = &fasta.Fasta{Name: "toad", Seq: dna.StringToBases("TTGTTATTC")}
var ahsoka = &fasta.Fasta{Name: "ahsoka", Seq: dna.StringToBases("TTGTTC")}

//make needed cigar and graph using the graph function being tested
var _, aln = align.ConstGap(toad.Seq, ahsoka.Seq, align.HumanChimpTwoScoreMatrix, -430)
var graph = cigarToGraph(toad, ahsoka, aln)

//make a manual graph using the fake fastas
var manual = buildManualGraph(toad, ahsoka)

//Did graph and manual end up with the same number of nodes?
func TestAddNode(t *testing.T) {
	if len(graph.Nodes) != len(aln) {
		t.Errorf("Number of nodes does not match with the cigar")
	}
}

//test if graph and manual are equivalent.
func TestCigarToGraphMatchManual(t *testing.T) {
	compareScore := 0
	seqCounter := compareSeq(manual, graph)
	edgeCounter := compareEdges(manual, graph)
	compareScore = seqCounter + edgeCounter

	if compareScore != 0 {
		t.Errorf("'Graph' was not equivalent to 'manual'. Expect counter value to be zero. Use location of non-zero output to determine the locatino of the issue.  compareSeq() counter was %v, compareEdge() counter was %v, giving a total score of %v", seqCounter, edgeCounter, compareScore)
	}
}

//functions below this line were used to execute test functions

//this function makes the manually created graph (called 'manual')
func buildManualGraph(toad *fasta.Fasta, ahsoka *fasta.Fasta) *simpleGraph.SimpleGraph {
	manual := simpleGraph.NewGraph()
	nodeOne := &simpleGraph.Node{
		Id:   0,
		Name: "toad_0_3",
		Seq:  toad.Seq[0:3]}
	nodeTwo := &simpleGraph.Node{
		Id:   1,
		Name: "toad_3_6",
		Seq:  toad.Seq[3:6]}
	nodeThree := &simpleGraph.Node{
		Id:   2,
		Name: "toad_6_9",
		Seq:  toad.Seq[6:9]}

	simpleGraph.AddNode(manual, nodeOne)
	simpleGraph.AddNode(manual, nodeTwo)
	simpleGraph.AddEdge(nodeOne, nodeTwo, 0.5)
	simpleGraph.AddNode(manual, nodeThree)
	simpleGraph.AddEdge(nodeTwo, nodeThree, 1.0)
	simpleGraph.AddEdge(nodeOne, nodeThree, 0.5)

	return manual
}

//this function compares the sequences at each node, going base by base.
func compareSeq(manual *simpleGraph.SimpleGraph, graph *simpleGraph.SimpleGraph) int {
	var seqCounter int = 0
	for i := 0; i < len(manual.Nodes); i++ {
		for j := 0; j < len(manual.Nodes[i].Seq); j++ {
			if manual.Nodes[i].Seq[j] == graph.Nodes[i].Seq[j] {
				seqCounter += 0
			} else {
				seqCounter += 1
			}
		}
	}
	return seqCounter
}

//this function compares the edges (both previous and next) by comparing the name of the node within the prev[] and next[] for each node. As a reminder, the prev[] and next[] are type *Node, so they have all the information of the connecting nodes.
//TODO: Compare the sequences within each edge between the 'graph' and 'manual'.
//TODO: Compare the probability of each edge between the 'graph' and 'manual'.
func compareEdges(manual *simpleGraph.SimpleGraph, graph *simpleGraph.SimpleGraph) int {
	totalManualPrevEdges := 0
	totalManualNextEdges := 0

	for i := 0; i < len(manual.Nodes); i++ {
		totalManualPrevEdges += len(manual.Nodes[i].Prev)
	}
	for i := 0; i < len(manual.Nodes); i++ {
		totalManualNextEdges += len(manual.Nodes[i].Next)
	}

	prevNodeMatchCounter := 0
	for i := 1; i < len(manual.Nodes); i++ {
		for j := 0; j < len(manual.Nodes[i].Prev); j++ {
			if len(manual.Nodes[i].Prev[j].Dest.Seq) == len(graph.Nodes[i].Prev[j].Dest.Seq) {
				prevNodeMatchCounter += 0
			} else {
				prevNodeMatchCounter += 1
			}
		}
	}

	nextNodeMatchCounter := 0
	for i := 0; i < len(manual.Nodes)-1; i++ {
		for j := 0; j < len(manual.Nodes[i].Next); j++ {
			if len(manual.Nodes[i].Next[j].Dest.Seq) == len(graph.Nodes[i].Next[j].Dest.Seq) {
				nextNodeMatchCounter += 0
			} else {
				nextNodeMatchCounter += 1
			}
		}
	}

	totalEdgeCounter := prevNodeMatchCounter + nextNodeMatchCounter

	return totalEdgeCounter
}
