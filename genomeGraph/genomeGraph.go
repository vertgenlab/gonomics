// Package genomeGraph has structs and tools for reading, writing, editing and aligning graph representations of genomes
package genomeGraph

import (
	"fmt"
	"io"
	"log"
	"math"
	"strings"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
)

// GenomeGraph struct contains a slice of Nodes.
type GenomeGraph struct {
	Nodes []Node
}

// Node is uniquely definded by Id and is encoded with information
// describing sequence order and orientation and annotated variance.
type Node struct {
	Id        uint32
	ColumnId  uint32
	Seq       []dna.Base        // only this field or the SeqThreeBit will be kept
	SeqTwoBit *dnaTwoBit.TwoBit // this will change to a ThreeBit or be removed
	Prev      []Edge
	Next      []Edge
} // used to have Name (string) and Info (Annotation)

// Edge describes the neighboring nodes and a weighted probability
// of the more likely path.
type Edge struct {
	Dest *Node
	Prob float32
}

// Annotation struct is an uint64 encoding of allele id, starting position on linear reference and variant on node
// a single byte will represent what allele the node came from, uint32 will be used for starting position of chromosome of the linear reference
// uint8 are variants are represented as follows: 0=match, 1=mismatch, 2=insert, 3=deletion, 4=hap
/*type Annotation struct {
	Start   uint32
	Allele  byte
	Variant uint8
}*/

// Read will process a simple graph formated text file and parse the data into graph fields.
func Read(filename string) *GenomeGraph {
	simpleReader := fileio.NewByteReader(filename)
	genome := EmptyGraph()
	var weight float32
	var line string
	var words []string = make([]string, 0, 2)
	var nodeId, homeNodeIdx, destNodeIdx uint32
	var homeNode, destNode *Node
	var i int

	for reader, done := fileio.ReadLine(simpleReader); !done; reader, done = fileio.ReadLine(simpleReader) {
		line = reader.String()
		switch true {
		case strings.HasPrefix(line, ">"):
			nodeId = parse.StringToUint32(line[1:])
			AddNode(genome, &Node{Id: nodeId})
		case strings.Contains(line, "\t"):
			words = strings.Split(line, "\t")
			homeNodeIdx = parse.StringToUint32(words[0])
			homeNode = &genome.Nodes[homeNodeIdx]
			if len(words) > 2 {
				for i = 1; i < len(words); i += 2 {
					weight = parse.StringToFloat32(words[i])
					destNodeIdx = parse.StringToUint32(words[i+1])
					destNode = &genome.Nodes[destNodeIdx]
					AddEdge(homeNode, destNode, weight)
				}
			}
		default:
			genome.Nodes[nodeId].Seq = append(genome.Nodes[nodeId].Seq, dna.ByteSliceToDnaBases(reader.Bytes())...)
		}
	}
	for i = range genome.Nodes {
		if len(genome.Nodes[i].Seq) != 0 {
			genome.Nodes[i].SeqTwoBit = dnaTwoBit.NewTwoBit(genome.Nodes[i].Seq)
		}
	}
	return genome
}

// AddNode will add the values in n to the graph at the index of n.Id
// A pointer to the new location of the node (inside the graph) is returned.
func AddNode(g *GenomeGraph, n *Node) *Node {
	const positionsToExtend = 1000 // when we need to increase the slice, do this many nodes at a time
	if int(n.Id) < len(g.Nodes) {  // if the memory for this node has already been allocated
		// then we can overwrite it as long as it does not already exist
		if len(g.Nodes[n.Id].Seq) != 0 { // if the node already exists because data has been written here
			log.Panicf("Error: tried to add a node of id=%d, when that id already exists\n", n.Id)
		}
	} else if int(n.Id) < cap(g.Nodes) { // if we already have the capacity, but not the length
		g.Nodes = g.Nodes[:n.Id+1]
	} else { // if we need to increase capacity
		futureNodes := make([]Node, int(n.Id)+1, numbers.Max(cap(g.Nodes)+positionsToExtend, int(n.Id)+1))
		copy(futureNodes, g.Nodes)
		g.Nodes = futureNodes
	}
	g.Nodes[n.Id] = *n
	return &g.Nodes[n.Id]
}

// AddEdge will append two edges one forward and one backwards for any two
// given node. Provide a probability float32 to specify a weight for an edge
// to describe the more likely path through the graph.
func AddEdge(u, v *Node, p float32) {
	u.Next = append(u.Next, Edge{Dest: v, Prob: p})
	v.Prev = append(v.Prev, Edge{Dest: u, Prob: p})
}

// SetEvenWeights will loop through a slice of edges and set the probability weight
// divided by the length of the slice.
func SetEvenWeights(u *Node) {
	var edge int
	var weights float32 = 1 / float32(len(u.Next))
	for edge = 0; edge < len(u.Next); edge++ {
		u.Next[edge].Prob = weights
	}
}

// Write function will process GenomeGraph and write the data to a file.
func Write(filename string, sg *GenomeGraph) {
	lineLength := 50
	file := fileio.EasyCreate(filename)
	defer file.Close()
	WriteToGraphHandle(file, sg, lineLength)
}

// EmptyGraph will allocate a new zero pointer to a simple graph and will allocate memory for the Nodes of the graph.
func EmptyGraph() *GenomeGraph {
	return &GenomeGraph{Nodes: make([]Node, 0)}
}

// PrintGraph will quickly print simpleGraph to standard out.
func PrintGraph(gg *GenomeGraph) {
	Write("/dev/stdout", gg)
}

// WriteToGraphHandle will help with any error handling when writing GenomeGraph to file.
func WriteToGraphHandle(file io.Writer, gg *GenomeGraph, lineLength int) {
	var err error
	var i, j int
	for i = 0; i < len(gg.Nodes); i++ {
		_, err = fmt.Fprintf(file, ">%d\n", gg.Nodes[i].Id)
		exception.PanicOnErr(err)
		for j = 0; j < len(gg.Nodes[i].Seq); j += lineLength {
			if j+lineLength > len(gg.Nodes[i].Seq) {
				_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(gg.Nodes[i].Seq[j:]))
				exception.PanicOnErr(err)
			} else {
				_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(gg.Nodes[i].Seq[j:j+lineLength]))
				exception.PanicOnErr(err)
			}
		}
	}
	for i = 0; i < len(gg.Nodes); i++ {
		near := gg.Nodes[i].Next
		if len(near) > 0 {
			_, err = fmt.Fprintf(file, "%d", gg.Nodes[i].Id)
			exception.PanicOnErr(err)
			for j = 0; j < len(near); j++ {
				_, err = fmt.Fprintf(file, "\t%v\t%d", near[j].Prob, near[j].Dest.Id)
				exception.PanicOnErr(err)
			}
			_, err = fmt.Fprintf(file, "\n")
			exception.PanicOnErr(err)
		}
	}
}

// BasesInGraph will calculate the number of bases contained in GenomeGraph using dnaTwoBit.
func BasesInGraph(g *GenomeGraph) int {
	var i, baseCount int = 0, 0
	for i = 0; i < len(g.Nodes); i++ {
		baseCount += len(g.Nodes[i].Seq)
	}
	return baseCount
}

// AlignGraphTraversal performs graph traversal for alignment either to the left or right, calculating optimal alignment using dynamic programming.
func AlignGraphTraversal(n *Node, seq []dna.Base, position int, currentPath []uint32, extension int, read []dna.Base, scores [][]int64, matrix *MatrixAln, sk *scoreKeeper, dynamicScore *dynamicScoreKeeper, pool *sync.Pool, direction byte) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	// Fetch a pooled object to reduce allocations.
	s := pool.Get().(*dnaPool)
	defer pool.Put(s)

	// Initialize sequence and path slices without reallocation.
	s.Seq = append(s.Seq[:0], seq...)
	s.Path = append(s.Path[:0], currentPath...)

	// Traverse node sequence and update the path accordingly.
	s.Seq = nodeSeqTraversal(n, extension, position, seq, s.Seq, direction)
	if direction == leftTraversal {
		s.Path = append(s.Path, n.Id) // Add node ID for left direction.
	}

	// Initialize local variables to hold alignment details.
	var alignment []cigar.ByteCigar
	var score int64
	var targetStart, queryStart int

	// Check termination condition for recursion.
	isEndOfTraversal := (direction == leftTraversal && (len(s.Seq) >= extension || len(n.Prev) == 0)) ||
		(direction == rightTraversal && (len(s.Seq) >= extension || len(n.Next) == 0))

	if isEndOfTraversal {
		// Perform dynamic alignment at the end of traversal.
		score, alignment, targetStart, queryStart = DynamicAln(s.Seq, read, scores, matrix, -600, dynamicScore, direction)
		if direction == rightTraversal {
			targetStart += position // Adjust target start for right direction.
		}
		return alignment, score, targetStart, queryStart, append([]uint32(nil), s.Path...)
	}

	// Recursive case: traverse next or previous nodes based on direction.
	score = math.MinInt64
	nodes := n.Next
	if direction == leftTraversal {
		nodes = n.Prev
	}
	for _, edge := range nodes {
		route, currScore, end, start, path := AlignGraphTraversal(edge.Dest, s.Seq, 0, s.Path, extension, read, scores, matrix, sk, dynamicScore, pool, direction)
		if currScore > score {
			// Update alignment details if current score is better.
			score, alignment, targetStart, queryStart = currScore, route, end, start
			s.Path = path
		}
	}

	// Adjust the target start for left direction after all recursive calls.
	if direction == leftTraversal {
		targetStart = position - len(s.Seq) + targetStart
		cigar.ReverseBytesCigar(alignment)
		ReversePath(s.Path)
	} else {
		targetStart += position // Adjust for right direction.
	}

	return alignment, score, targetStart, queryStart, append([]uint32(nil), s.Path...)
}

// DynamicAln performs dynamic programming-based sequence alignment by traversing a graph of nodes
func DynamicAln(alpha, beta []dna.Base, scores [][]int64, matrix *MatrixAln, gapPen int64, dynamicScore *dynamicScoreKeeper, direction byte) (int64, []cigar.ByteCigar, int, int) {
	dynamicScore.currMax = 0 // Reset the current max score.
	var i, j int
	var rows, columns int = len(alpha), len(beta)
	// Initialize the first row and column based on gap penalties.
	for i = 1; i <= rows; i++ {
		matrix.m[i][0] = int64(i) * gapPen
		matrix.trace[i][0] = 'D' // Indicate a deletion.
	}
	for j = 1; j <= columns; j++ {
		matrix.m[0][j] = int64(j) * gapPen
		matrix.trace[0][j] = 'I' // Indicate an insertion.
	}

	// Fill in the scoring and traceback matrices.
	for i = 1; i <= rows; i++ {
		for j = 1; j <= columns; j++ {
			matchOrMismatch := matrix.m[i-1][j-1] + scores[alpha[i-1]][beta[j-1]]
			deletion := matrix.m[i-1][j] + gapPen
			insertion := matrix.m[i][j-1] + gapPen

			// Choose the best score.
			if matchOrMismatch >= deletion && matchOrMismatch >= insertion {
				matrix.m[i][j] = matchOrMismatch
				matrix.trace[i][j] = cigar.Match
			} else if deletion > insertion {
				matrix.m[i][j] = deletion
				matrix.trace[i][j] = cigar.Deletion
			} else {
				matrix.m[i][j] = insertion
				matrix.trace[i][j] = cigar.Insertion
			}
			// Update the current maximum score.
			if matrix.m[i][j] > dynamicScore.currMax {
				dynamicScore.currMax = matrix.m[i][j]
				dynamicScore.i, dynamicScore.j = i, j
			}
		}
	}

	// Traceback from the maximum score's position to build the alignment.
	alignment := make([]cigar.ByteCigar, 0)
	i, j = dynamicScore.i, dynamicScore.j
	for i > 0 && j > 0 {
		switch matrix.trace[i][j] {
		case cigar.Match:
			alignment = cigar.AddCigarByte(alignment, cigar.ByteCigar{RunLen: 1, Op: cigar.Match})
			i--
			j--
		case cigar.Insertion:
			alignment = cigar.AddCigarByte(alignment, cigar.ByteCigar{RunLen: 1, Op: cigar.Insertion})
			j--
		case cigar.Deletion:
			alignment = cigar.AddCigarByte(alignment, cigar.ByteCigar{RunLen: 1, Op: cigar.Deletion})
			i--
		}
	}
	// Prepend soft clips if necessary (not shown here, depends on alignment strategy).
	return dynamicScore.currMax, alignment, dynamicScore.i, dynamicScore.j
}

func nodeSeqTraversal(n *Node, extension int, position int, seq, ans []dna.Base, direction byte) []dna.Base {
	switch direction {
	case leftTraversal:
		startPos := position - extension
		if startPos < 0 {
			startPos = 0
		}
		// Ensure the combined sequence does not exceed the original sequence length.
		endPos := position
		if endPos > len(n.Seq) {
			endPos = len(n.Seq)
		}
		// Combine sequences.
		ans = append(ans, n.Seq[startPos:endPos]...)
		ans = append(ans, seq...)
	case rightTraversal:
		endPos := position + extension
		if endPos > len(n.Seq) {
			endPos = len(n.Seq)
		}
		// Combine sequences.
		ans = append(ans, seq...)
		ans = append(ans, n.Seq[position:endPos]...)
	}
	return ans
}
