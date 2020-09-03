package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

// SimpleGraph struct contains a slice of Nodes
type SimpleGraph struct {
	Nodes []*Node
}

// Node is uniquely definded by Id and is encoded with information
// describing sequence order and orientation and annotated variance
type Node struct {
	Id        uint32
	Name      string
	Seq       []dna.Base
	SeqTwoBit *dnaTwoBit.TwoBit
	Prev      []*Edge
	Next      []*Edge
	Info      Annotation
}

// Edge describes the neighboring nodes and a weighted probabilty
// of the of the more likely path
type Edge struct {
	Dest *Node
	Prob float32
}

// Annotation struct is an uint64 encoding of allele id, starting position on linear reference and variant on node
// a single byte will represent what allele the node came from, uint32 will be used for starting postion of chromosome of the linear reference
// uint8 are variants are represented as follows: 0=match, 1=mismatch, 2=insert, 3=deletion, 4=hap
type Annotation struct {
	Start   uint32
	Allele  byte
	Variant uint8
}

// Read will process a simple graph formated text file and parse the data into graph fields
func Read(filename string) *SimpleGraph {
	simpleReader := fileio.NewSimpleReader(filename)
	genome := NewGraph()
	var currNode *Node
	var edges map[string]*Node = make(map[string]*Node)
	var weight float32
	var line string
	var words []string = make([]string, 0, 2)
	var text []string = make([]string, 0, 3)
	var seqIdx int32 = -1
	var i int
	var ok bool

	for reader, done := fileio.ReadLine(simpleReader); !done; reader, done = fileio.ReadLine(simpleReader) {
		line = reader.String()
		switch true {
		case strings.HasPrefix(line, ">"):
			seqIdx++
			words = strings.Split(line[1:], ":")

			currNode = &Node{Id: uint32(seqIdx), Name: string(words[0])}
			if len(words) == 2 {
				text = strings.Split(words[1], "_")
				currNode.Info = Annotation{Allele: uint8(common.StringToUint32(text[1])), Start: common.StringToUint32(text[3]), Variant: uint8(common.StringToUint32(text[2]))}
			}
			AddNode(genome, currNode)
			_, ok = edges[line[1:]]
			if !ok {
				edges[string(line[1:])] = currNode
			}
		case strings.Contains(line, "\t"):
			words = strings.Split(line, "\t")
			if len(words) > 2 {
				for i = 1; i < len(words); i += 2 {
					weight = float32(common.StringToFloat64(words[i]))
					AddEdge(edges[words[0]], edges[words[i+1]], weight)
				}
			}
		case !strings.ContainsAny(line, "\t:"):

			genome.Nodes[seqIdx].Seq = append(genome.Nodes[seqIdx].Seq, dna.ByteSliceToDnaBases(reader.Bytes())...)
		}
	}
	for i = 0; i < len(genome.Nodes); i++ {
		genome.Nodes[i].SeqTwoBit = dnaTwoBit.NewTwoBit(genome.Nodes[i].Seq)
	}
	return genome
}

// ReadToMap will process a text file and create a hash look up using
// name string as keys to return a SimpleGraph. One possible use case
// is using this function to split the graph by chromocomes and using
// the name of the chromosomes as keys.
func ReadToMap(filename string) map[string]*SimpleGraph {
	genomeGraph := make(map[string]*SimpleGraph)
	var chrGraph string
	var line string
	var currSeq []dna.Base
	//var nodeId int64 = -1
	var nodeId uint32
	var doneReading bool = false
	var words, text []string
	var weight float32
	//var cInfo chromInfo.ChromInfo
	file := fileio.EasyOpen(filename)
	defer file.Close()
	//creates map: name points to Node
	//uses this map to add edges to graph
	edges := make(map[string]*Node)
	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, ">") {
			words = strings.Split(line, ":")
			chrGraph = words[0][1:]
			_, key := genomeGraph[chrGraph]
			if !key {
				genomeGraph[chrGraph] = NewGraph()
				nodeId = 0
			} else {
				nodeId++
			}
			tmp := Node{Id: nodeId, Name: words[0][1:]}
			if len(words) == 2 {
				text = strings.Split(words[1], "_")
				tmp.Info = Annotation{Allele: uint8(common.StringToUint32(text[1])), Start: common.StringToUint32(text[3]), Variant: uint8(common.StringToUint32(text[2]))}
			}
			AddNode(genomeGraph[chrGraph], &tmp)
			_, ok := edges[line[1:]]
			if !ok {
				edges[line[1:]] = &tmp
			}
		} else if strings.Contains(line, "\t") {
			words = strings.Split(line, "\t")
			if len(words) > 2 {
				for i := 1; i < len(words); i += 2 {
					weight = float32(common.StringToFloat64(words[i]))
					AddEdge(edges[words[0]], edges[words[i+1]], weight)
				}
			}
		} else {
			if !strings.Contains(line, ":") && !strings.Contains(line, "\t") {
				currSeq = dna.StringToBases(line)
				dna.AllToUpper(currSeq)
				genomeGraph[chrGraph].Nodes[nodeId].Seq = append(genomeGraph[chrGraph].Nodes[nodeId].Seq, currSeq...)
			}
		}
	}
	return genomeGraph
}

// AddNode will append a new Node to a slice of nodes in SimpleGraph
func AddNode(g *SimpleGraph, n *Node) {
	g.Nodes = append(g.Nodes, n)
}

// AddEdge will append two edges one forward and one backwards for any two
// given node. Provide a probability float32 to specify a weight for an edge
// to describe the more likely path through the graph
func AddEdge(u, v *Node, p float32) {
	u.Next = append(u.Next, &Edge{Dest: v, Prob: p})
	v.Prev = append(v.Prev, &Edge{Dest: u, Prob: p})
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

// Write function will process SimpleGraph and write the data to a file
func Write(filename string, sg *SimpleGraph) {
	lineLength := 50
	file := fileio.EasyCreate(filename)
	defer file.Close()
	WriteToGraphHandle(file, sg, lineLength)
}

// NewGraph will allocate a new zero pointer to a simple graph and will allocate memory for the Nodes of the graph
func NewGraph() *SimpleGraph {
	graph := new(SimpleGraph)
	graph.Nodes = make([]*Node, 0)
	return graph
}

// PrintGraph will quickly print simpleGraph to standard out
func PrintGraph(gg *SimpleGraph) {
	Write("/dev/stdout", gg)
}

// WriteToGraphHandle will help with any error handling when writing SimpleGraph to file.
func WriteToGraphHandle(file io.Writer, gg *SimpleGraph, lineLength int) {
	var err error
	var i, j int
	for i = 0; i < len(gg.Nodes); i++ {
		if &gg.Nodes[i].Info != nil {
			_, err = fmt.Fprintf(file, ">%s:%d_%d_%d_%d\n", gg.Nodes[i].Name, gg.Nodes[i].Id, gg.Nodes[i].Info.Allele, gg.Nodes[i].Info.Variant, gg.Nodes[i].Info.Start)
		} else {
			_, err = fmt.Fprintf(file, ">%s\n", gg.Nodes[i].Name)
		}
		for j = 0; j < len(gg.Nodes[i].Seq); j += lineLength {
			if j+lineLength > len(gg.Nodes[i].Seq) {
				_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(gg.Nodes[i].Seq[j:]))
				common.ExitIfError(err)
			} else {
				_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(gg.Nodes[i].Seq[j:j+lineLength]))
				common.ExitIfError(err)
			}
		}
	}
	for i = 0; i < len(gg.Nodes); i++ {
		if &gg.Nodes[i].Info != nil {
			_, err = fmt.Fprintf(file, "%s:%d_%d_%d_%d", gg.Nodes[i].Name, gg.Nodes[i].Id, gg.Nodes[i].Info.Allele, gg.Nodes[i].Info.Variant, gg.Nodes[i].Info.Start)
		} else {
			_, err = fmt.Fprintf(file, "%s", gg.Nodes[i].Name)
		}
		near := gg.Nodes[i].Next
		for j = 0; j < len(near); j++ {
			if &gg.Nodes[i].Info != nil {
				_, err = fmt.Fprintf(file, "\t%v\t%s:%d_%d_%d_%d", near[j].Prob, near[j].Dest.Name, near[j].Dest.Id, near[j].Dest.Info.Allele, near[j].Dest.Info.Variant, near[j].Dest.Info.Start)
			} else {
				_, err = fmt.Fprintf(file, "\t%v\t%s", near[j].Prob, near[j].Dest.Name)
			}
			common.ExitIfError(err)
		}
		_, err = fmt.Fprintf(file, "\n")
		common.ExitIfError(err)
	}
}

// WriteToGraphSplit will write all SimpleGraphs contained in a provided hash map to separate
// files based on unique keys.
func WriteToGraphSplit(filename string, gg map[string]*SimpleGraph) {
	var name string
	for chr := range gg {
		name = chr + "_" + filename + ".gg"
		Write(name, gg[chr])
	}
}

// BasesInGraph will calculate the number of bases contained in SimpleGraph using dnaTwoBit
func BasesInGraph(g *SimpleGraph) int {
	var i, baseCount int = 0, 0
	for i = 0; i < len(g.Nodes); i++ {
		baseCount += g.Nodes[i].SeqTwoBit.Len
	}
	return baseCount
}

// GraphToMap converts a genome graph to a map leading to Nodes keyed by the node names
func GraphToMap(g *SimpleGraph) map[string]*Node {
	m := make(map[string]*Node)
	for _, val := range g.Nodes {
		if _, ok := m[val.Name]; ok {
			log.Fatalf("ERROR: Node names in graph are not unique: %s", val.Name)
		} else {
			m[val.Name] = val
		}
	}
	return m
}
