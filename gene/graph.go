package gene

import (
	"errors"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/genomeGraph"
)

// GeneGraph defines a path through GeneGraph.Graph that corresponds to the linear sequence
// in the embedded Gene. The GeneGraph.Path is represented as a new GenomeGraph that copies the
// sequence in Graph for all nodes encompassing the sequence in Gene. Each node in Path has a
// maximum of 1 Next or Prev nodes.
//
// Nodes in Path are an identical copy of a node in Graph except
// the Node.Seq field may contain a subset of the sequence in the original node if the beginning
// or end of the Gene sequence does not align with node boundaries
// (i.e. Path.Node.Seq = Graph.Node.Seq[Gene.Start : Gene.End])
//
// len(Path.Nodes) is equal to len(Graph.Nodes), however indexes in Path.Nodes may be nil
// if the corresponding node in Graph.Nodes does not contain any sequence from Gene.
//
// PathNodeIds contains the Id for each non-nil node present in Path in order of occurrence
// in the Gene sequence. (i.e. Path.Nodes[PathNodeIds[0]].Next == Path.Nodes[PathNodeIds[1]])
//
// Note that all non-nil Nodes in Path (including Nodes with a subset of Seq) share an underlying
// Seq array with the corresponding Node in Graph. Care should be taken when trying to manipulate
// the Seq of a Node in Path.
type GeneGraph struct {
	*Gene
	PathNodeIds []uint32                // Id for each non-nil node in Path
	Path        genomeGraph.GenomeGraph // Sub-GenomeGraph where each node has max 1 Next/Prev
	Graph       genomeGraph.GenomeGraph // Complete GenomeGraph with all Next and Prev nodes
}

var MatchingPathNotFound = errors.New("error: matching path not found")

// SeqToPath finds a path through a graph that corresponds to the input
// sequence starting from startNode.Seq[startNodePos] and following edges.
// Sequence must be a perfect match (case ignored) to qualify as a valid path.
//
// Path is reported as a new GenomeGraph where each node has exactly 1 next/prev
// (expect for first and last node). error != nil if and only if a path through
// the graph with an identical sequence could not be found. On a non-nil error return
// the path up to the node where a mismatch occurred will be returned.
//
// NOTE: The first and last nodes in the path may contain a subset of the
// seq field (no change to seqTwoBit) if the input seq does not encompass the
// entire node sequence such that Node[0].Seq[0] is the first base of the input
// sequence and Node[last].Seq[last] is the last base of the input sequence.
func SeqToPath(seq []dna.Base, startNode *genomeGraph.Node, startNodePos int, graph genomeGraph.GenomeGraph) (ids []uint32, answer genomeGraph.GenomeGraph, err error) {
	answer = genomeGraph.GenomeGraph{Nodes: make([]*genomeGraph.Node, len(graph.Nodes))}

	// NOTE: currNode is always a copy of the original node with a subset path
	// for readability, all calls to a node from the full graph are done either
	// through graph.Nodes[i] or the variable 'origNode'
	var currNode *genomeGraph.Node = copyNodeWithSubPath(startNode, nil) // start graph with first node
	currNode.Seq = currNode.Seq[startNodePos:]                           // subset seq node in copy

	currSeqIdx := len(startNode.Seq[startNodePos:])

	// check startNode is valid
	if dna.CompareSeqsIgnoreCase(seq[:currSeqIdx], startNode.Seq[startNodePos:]) != 0 {
		err = MatchingPathNotFound
		return
	} else { // add startNode copy to answer
		answer.Nodes[currNode.Id] = currNode
		ids = append(ids, currNode.Id)
	}

	var endIdx int // this is the index of the last base in the nodes seq that is consumed by input seq
	for ; currSeqIdx < len(seq); currSeqIdx += len(currNode.Seq) {

		currNode, endIdx, err = findNodeInNext(seq[currSeqIdx:], graph.Nodes[currNode.Id], currNode)
		if err != nil {
			return
		}
		currNode.Seq = currNode.Seq[:endIdx] // only changes seq of last node if not fully consumed by input seq

		answer.Nodes[currNode.Id] = currNode // set currNode in graph
		ids = append(ids, currNode.Id)
	}
	return
}

// findNodeInNext finds the next node with an identical sequence to the input seq.
func findNodeInNext(seq []dna.Base, origNode *genomeGraph.Node, currNode *genomeGraph.Node) (nextNode *genomeGraph.Node, endIdx int, err error) {
	var seqToTest []dna.Base
	for _, e := range origNode.Next {
		if len(e.Dest.Seq) > len(seq) {
			seqToTest = seq
		} else {
			seqToTest = seq[:len(e.Dest.Seq)]
		}

		if dna.CompareSeqsIgnoreCase(seqToTest, e.Dest.Seq[:len(seqToTest)]) == 0 {
			nextNode = copyNodeWithSubPath(e.Dest, currNode)
			endIdx = len(seqToTest)
			return
		}
	}
	err = MatchingPathNotFound
	return
}

// copyNodeWithSubPath will copy the input node and and set next and prev.
func copyNodeWithSubPath(origNode *genomeGraph.Node, currNode *genomeGraph.Node) *genomeGraph.Node {
	nextNode := new(genomeGraph.Node)
	*nextNode = *origNode // copy all fields in input origNode. note that slices are referenced

	// remove original Next and Prev values
	nextNode.Prev = nil
	nextNode.Next = nil

	if currNode != nil { // set curr.Next and next.Prev
		currNode.Next = []genomeGraph.Edge{{Dest: nextNode, Prob: 1}}
		nextNode.Prev = []genomeGraph.Edge{{Dest: currNode, Prob: 1}}
	}

	return nextNode
}
