package graphReconstruct

import (
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/genomeGraph"
)

type graphColumn struct {
	AlignId    int
	AlignNodes map[string][]*genomeGraph.Node //string keys refer to species that key to a slice of pointers to the nodes of that species that fall into a single slignment column
}

// BuildNodes uses a graphColumn to create nodes for an ancestor's graph seq that represents all the unique sequences in an aligned graph.
func BuildNodes(root *expandedTree.ETree, column graphColumn, id uint32) uint32 {
	var nodeInfo = make(map[string]bool)
	for _, nodes := range column.AlignNodes { //nodes is all nodes for an individual species
		for n := range nodes { //n is an individual node of an individual species
			stringSeq := dna.BasesToString(nodes[n].Seq)
			nodeInfo[stringSeq] = true
		}
	}
	for seq := range nodeInfo {
		var newNode *genomeGraph.Node
		newNode = &genomeGraph.Node{Id: id, Seq: dna.StringToBases(seq), SeqTwoBit: dnaTwoBit.NewTwoBit(dna.StringToBases(seq)), Next: nil, Prev: nil}
		column.AlignNodes[root.Name] = append(column.AlignNodes[root.Name], newNode)
		id += 1
	}
	return id
}

//BuildEdges connects the nodes of a species' graph that are stored in GraphColumns
//func BuildEdges
//start without prob
//loop through species in column, go through all nodes
//check if that node's seq matches the seq of the ancestor node without an edge
//make the next of that node the same as this node's next
//FindAncSeq creates a graph from the node records stored in GraphColumns and then calls PathFinder and seqOfPath to determine the most likley seq of the ancestor before assigning that
//seq to the Fasta field of the ancestors tree node
//func FindAncSeq will loop through aligncolumns and build a single graph of all of the nodes that belong to the ancestor species after edges are created
//run PathFinder on the graph for the anc, run seqOfPath, then turn that to a fasta for that node of the tree

// seqOfPath takes in a graph and a path specified by the Node IDs and returns the seq of the path through the graph.
func seqOfPath(g *genomeGraph.GenomeGraph, path []uint32) []dna.Base {
	var seq []dna.Base
	var foundInGraph = false
	for p := 0; p < len(path); p++ {
		foundInGraph = false
		for n := 0; n < len(g.Nodes) && !foundInGraph; n++ {
			if g.Nodes[n].Id == path[p] {
				foundInGraph = true
				seq = append(seq, g.Nodes[n].Seq...)
			} else {
			}
		}
		if !foundInGraph {
			log.Fatal("path is invalid")
		}
	}
	return seq
}

// PathFinder takes a graph and returns the most likely path through that graph after checking all possible paths from the first node to the last.
func PathFinder(g *genomeGraph.GenomeGraph) ([]uint32, float32) {
	var finalPath []uint32
	var finalProb float32
	var tempPath = make([]uint32, 0)

	for n := 0; n < len(g.Nodes); n++ {
		if g.Nodes[n].Id == 0 {
			finalProb, finalPath = bestPath(&g.Nodes[n], 1, tempPath)
		}
	}
	return finalPath, finalProb
}

// bestPath is the helper function for PathFinder, and recursively traverses the graph depth first to determine the most likely path from start to finish.
func bestPath(node *genomeGraph.Node, prevProb float32, path []uint32) (prob float32, pathOut []uint32) {
	var tempProb float32 = 0
	var finalProb float32
	var finalPath []uint32

	path = append(path, node.Id)
	if len(node.Next) == 0 {
		return prevProb, path
	}
	for i := range node.Next {
		tempProb = node.Next[i].Prob * prevProb
		currentProb, currentPath := bestPath(node.Next[i].Dest, tempProb, path)
		if currentProb > finalProb {
			finalProb = currentProb
			finalPath = currentPath
		}
	}
	return finalProb, finalPath
}
