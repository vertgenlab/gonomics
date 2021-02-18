package graph

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
)

// Graph is a set of nodes that together form a directed acyclic graph.
// Not all nodes are necessarily connected (i.e. graph may be discontiguous).
// Multiple graphs can exist with the same underlying data, but a distinct set of
// edges (e.g. distinct Human and Chimp Edge sets that point to the same
// sequence data, but may diverge and take distinct paths through the same
// underlying graph.
type Graph struct {
	Nodes []*Node // All Nodes present in graph. Ordered by Node.Idx
}

// Node is uniquely defined by Idx and contains a set of Edges defining possible
// paths through the graph.
type Node struct {
	Idx       uint32    // Index in Graph.Nodes. May be different for subgraphs
	Seq       *Data     // Sequence and associated identifiers
	Next      []Edge    // Possible paths forward in graph
	NextProbs []float32 // Probability of each Edge in Next. May be different for subgraphs
	Prev      []Edge    // Possible paths backwards in graph
	PrevProbs []float32 // Probability of each Edge in Prev. May be different for subgraphs
}

// Edge defines the path to another Node in the graph
type Edge *Node

// Data contains sequence order and orientation information as well as annotated variance
type Data struct {
	Name      string // Arbitrary name for data in struct (e.g. "LINE Element", "rs1011", etc). Not required to be unique.
	Seq       []dna.Base
	SeqTwoBit *dnaTwoBit.TwoBit
	Info      Annotation
}

// Annotation struct is an uint64 encoding of allele id, starting position on linear reference and variant on node
// a single byte will represent what allele the node came from, uint32 will be used for starting postion of chromosome of the linear reference
// uint8 are variants are represented as follows: 0=match, 1=mismatch, 2=insert, 3=deletion, 4=hap
type Annotation struct {
	Start   uint32
	Allele  byte
	Variant uint8
}
