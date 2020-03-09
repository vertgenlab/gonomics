package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"strings"
)

type SimpleGraph struct {
	Nodes []*Node
}

type Node struct {
	Id   uint32
	Name string
	Seq  []dna.Base
	Prev []*Edge
	Next []*Edge
	Info *Annotation
}

type Edge struct {
	Dest *Node
	Prob float32
}
//Storing annotation as a 64 uint: first 8 will be used to represent what allele, next 32 will be used for starting postion of chromosome.
//variants are represented as follows: 0=match, 1=mismatch, 2=insert, 3=deletion, 4=hap
type Annotation struct {
	Allele  uint8
	Start   uint32
	Variant uint8
}

func Read(filename string) *SimpleGraph {
	genomeGraph := NewGraph()
	var line string
	var currSeq []dna.Base
	var seqIdx int64 = -1
	var doneReading bool = false
	var words []string
	var weight float32
	file := fileio.EasyOpen(filename)
	defer file.Close()
	//creates map: name points to Node
	//uses this map to add edges to graph
	edges := make(map[string]*Node)
	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, ">") {
			seqIdx++
			words = strings.Split(line, "_")
			tmp := Node{Id: uint32(seqIdx), Name: words[0][1:]}
			if len(words) == 3 {
				tmp.Info = &Annotation{Start: common.StringToUint32(words[2])}
			}
			AddNode(genomeGraph, &tmp)
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
			if !strings.Contains(line, "_") {
				currSeq = dna.StringToBases(line)
				dna.AllToUpper(currSeq)
				genomeGraph.Nodes[seqIdx].Seq = append(genomeGraph.Nodes[seqIdx].Seq, currSeq...)
			}
		}
	}
	return genomeGraph
}

//TODO: We are transitioning to a new Read function that will keep track of chromosome lengths
func DevRead(filename string) (*SimpleGraph, map[string]*chromInfo.ChromInfo) {
	genomeGraph := NewGraph()
	chrSize := make(map[string]*chromInfo.ChromInfo)
	var count int64 = 0
	var line string
	var currSeq []dna.Base
	var seqIdx int64 = -1
	var doneReading bool = false
	var words []string
	var weight float32
	file := fileio.EasyOpen(filename)
	defer file.Close()
	//creates map: name points to Node
	//uses this map to add edges to graph
	edges := make(map[string]*Node)
	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, ">") {
			seqIdx++
			words = strings.Split(line, "_")
			tmp := Node{Id: uint32(seqIdx), Name: words[0][1:]}
			if len(words) == 3 {
				tmp.Info = &Annotation{Start: common.StringToUint32(words[2])}
			}
			AddNode(genomeGraph, &tmp)
			_, info := chrSize[words[0][1:]]
			if !info {
				cInfo := chromInfo.ChromInfo{Name: words[0][1:], Size: 0, Order: count}
				chrSize[words[0][1:]] = &cInfo
				count++
			}

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
			if !strings.Contains(line, "_") {
				currSeq = dna.StringToBases(line)
				dna.AllToUpper(currSeq)

				chrSize[genomeGraph.Nodes[seqIdx].Name].Size += int64(len(currSeq))
				genomeGraph.Nodes[seqIdx].Seq = append(genomeGraph.Nodes[seqIdx].Seq, currSeq...)
			}
		}
	}
	return genomeGraph, chrSize
}

func AddNode(g *SimpleGraph, n *Node) {
	g.Nodes = append(g.Nodes, n)
}

func AddEdge(u, v *Node, p float32) {
	u.Next = append(u.Next, &Edge{Dest: v, Prob: p})
	v.Prev = append(v.Prev, &Edge{Dest: u, Prob: p})
}

func SetEvenWeights(u *Node) {
	var edge int
	var weights float32 = 1 / float32(len(u.Next))
	for edge = 0; edge < len(u.Next); edge++ {
		u.Next[edge].Prob = weights
	}
	//do we care about the prev edges?
}

func Write(filename string, sg *SimpleGraph) {
	lineLength := 50
	file := fileio.EasyCreate(filename)
	defer file.Close()
	WriteToGraphHandle(file, sg, lineLength)
}

func NewGraph() *SimpleGraph {
	graph := new(SimpleGraph)
	graph.Nodes = make([]*Node, 0)
	return graph
}

func PrintGraph(gg *SimpleGraph) {
	Write("/dev/stdout", gg)
}

func WriteToGraphHandle(file io.Writer, gg *SimpleGraph, lineLength int) {
	var err error
	var i, j int
	for i = 0; i < len(gg.Nodes); i++ {
		if gg.Nodes[i].Info != nil {
			_, err = fmt.Fprintf(file, ">%s_%d_%d\n", gg.Nodes[i].Name, gg.Nodes[i].Id, gg.Nodes[i].Info.Start)
		} else {
			_, err = fmt.Fprintf(file, ">%s_%d\n", gg.Nodes[i].Name, gg.Nodes[i].Id)
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
		_, err = fmt.Fprintf(file, "%s_%d_%d", gg.Nodes[i].Name, gg.Nodes[i].Id, gg.Nodes[i].Info.Start)
		near := gg.Nodes[i].Next
		for j = 0; j < len(near); j++ {
			_, err = fmt.Fprintf(file, "\t%v\t%s_%d_%d", near[j].Prob, near[j].Dest.Name, near[j].Dest.Id, near[j].Dest.Info.Start)
			common.ExitIfError(err)
		}
		_, err = fmt.Fprintf(file, "\n")
		common.ExitIfError(err)
	}
}

func WriteToGraphCoordinates(file io.Writer, gg *SimpleGraph, lineLength int) {
	var err error
	var i, j int
	for i = 0; i < len(gg.Nodes); i++ {
		_, err = fmt.Fprintf(file, "%s\n", ">"+gg.Nodes[i].Name)
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
		_, err = fmt.Fprintf(file, "%s", gg.Nodes[i].Name)
		near := gg.Nodes[i].Next
		for j = 0; j < len(near); j++ {
			_, err = fmt.Fprintf(file, "\t%v\t%s", near[j].Prob, near[j].Dest.Name)
			common.ExitIfError(err)
		}
		_, err = fmt.Fprintf(file, "\n")
		common.ExitIfError(err)
	}
}
