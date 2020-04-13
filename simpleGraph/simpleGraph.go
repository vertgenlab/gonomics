package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"strings"
)

type SimpleGraph struct {
	Nodes []*Node
}

type Node struct {
	Id        uint32
	Name      string
	Seq       []dna.Base
	SeqTwoBit *dnaTwoBit.TwoBit
	Prev      []*Edge
	Next      []*Edge
	Info      *Annotation
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

//TODO: We are transitioning to a new Read function that will keep track of chromosome lengths
func Read(filename string) *SimpleGraph {
	genomeGraph := NewGraph()
	//chrSize := make(map[string]*chromInfo.ChromInfo)
	//var count int64 = 0
	var line string
	var currSeq []dna.Base
	var seqIdx int64 = -1
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
			seqIdx++
			words = strings.Split(line, ":")
			tmp := Node{Id: uint32(seqIdx), Name: words[0][1:], Info: nil}
			//_, cinfo := chrSize[tmp.Name]
			//if !cinfo {
			//	cInfo := chromInfo.ChromInfo{Name: tmp.Name, Size: 0, Order: count}
			//	chrSize[tmp.Name] = &cInfo
			//	count++
			//}
			if len(words) == 2 {
				text = strings.Split(words[1], "_")
				tmp.Info = &Annotation{Allele: uint8(common.StringToUint32(text[1])), Start: common.StringToUint32(text[3]), Variant: uint8(common.StringToUint32(text[2]))}
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
			if !strings.Contains(line, ":") && !strings.Contains(line, "\t") {
				currSeq = dna.StringToBases(line)
				dna.AllToUpper(currSeq)
				//if genomeGraph.Nodes[seqIdx].Info == nil {
				//	chrSize[genomeGraph.Nodes[seqIdx].Name].Size += int64(len(currSeq))
				//} else { //if genomeGraph.Nodes[seqIdx].Info != nil
				//	if genomeGraph.Nodes[seqIdx].Info.Allele == 0 {
				//		words = strings.Split(genomeGraph.Nodes[seqIdx].Name, ":")
				//		chrSize[words[0]].Size += int64(len(currSeq))
				//	}
				//}
				genomeGraph.Nodes[seqIdx].Seq = append(genomeGraph.Nodes[seqIdx].Seq, currSeq...)
			}
		}
	}
	for i := 0; i < len(genomeGraph.Nodes); i++ {
		genomeGraph.Nodes[i].SeqTwoBit = dnaTwoBit.NewTwoBit(genomeGraph.Nodes[i].Seq)
	}
	return genomeGraph
}

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
			tmp := Node{Id: nodeId, Name: words[0][1:], Info: nil}
			if len(words) == 2 {
				text = strings.Split(words[1], "_")
				tmp.Info = &Annotation{Allele: uint8(common.StringToUint32(text[1])), Start: common.StringToUint32(text[3]), Variant: uint8(common.StringToUint32(text[2]))}
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
		if gg.Nodes[i].Info != nil {
			_, err = fmt.Fprintf(file, "%s:%d_%d_%d_%d", gg.Nodes[i].Name, gg.Nodes[i].Id, gg.Nodes[i].Info.Allele, gg.Nodes[i].Info.Variant, gg.Nodes[i].Info.Start)
		} else {
			_, err = fmt.Fprintf(file, "%s", gg.Nodes[i].Name)
		}
		near := gg.Nodes[i].Next
		for j = 0; j < len(near); j++ {
			if gg.Nodes[i].Info != nil {
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

func WriteToGraphSplit(filename string, gg map[string]*SimpleGraph) {
	var name string
	for chr := range gg {
		name = chr + "_" + filename + ".gg"
		Write(name, gg[chr])
	}
}

func BasesInGraph(g *SimpleGraph) int {
	var i, baseCount int = 0, 0
	for i = 0; i < len(g.Nodes); i++ {
		baseCount += g.Nodes[i].SeqTwoBit.Len
	}
	return baseCount
}
