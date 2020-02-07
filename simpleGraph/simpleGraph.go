package simpleGraph

import (
	"fmt"
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
}

type Edge struct {
	Dest *Node
	Prob float64
}

func AddNode(g *SimpleGraph, n *Node) {
	g.Nodes = append(g.Nodes, n)
}

func AddEdge(u, v *Node, p float64) {
	u.Next = append(u.Next, &Edge{Dest: v, Prob: p})
	v.Prev = append(v.Prev, &Edge{Dest: u, Prob: p})
}

func Read(filename string) *SimpleGraph {
	answer := NewGraph()
	var line string
	//var name []string
	var currSeq []dna.Base
	var seqIdx int64 = -1
	var doneReading bool = false
	var words []string
	var weight float64
	file := fileio.EasyOpen(filename)
	defer file.Close()

	//creates map: name points to Node
	//uses this map to add edges to graph
	edges := make(map[string]*Node)

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, ">") {
			seqIdx++
			//name = strings.Split(line, ":")
			tmp := Node{Id: uint32(seqIdx), Name: line[1:]}
			AddNode(answer, &tmp)
			_, ok := edges[line[1:]]
			if !ok {
				edges[line[1:]] = &tmp
			}
		} else {
			if strings.Contains(line, "\t") {
				words = strings.Split(line, "\t")
				if len(words) > 2 {
					for i := 1; i < len(words); i += 2 {
						weight = common.StringToFloat64(words[i])
						AddEdge(edges[words[0]], edges[words[i+1]], weight)
					}
				}
			} else {
				if !strings.Contains(line, "_") {
					currSeq = dna.StringToBases(line)
					dna.AllToUpper(currSeq)
					answer.Nodes[seqIdx].Seq = append(answer.Nodes[seqIdx].Seq, currSeq...)
				}

			}
		}
	}
	return answer
}

func NewGraph() *SimpleGraph {
	graph := new(SimpleGraph)
	graph.Nodes = make([]*Node, 0)
	return graph
}

func WriteToFileHandle(file io.Writer, records []*Node, lineLength int) {
	var err error
	for _, rec := range records {
		_, err = fmt.Fprintf(file, ">%d\n", rec.Id)
		common.ExitIfError(err)
		for i := 0; i < len(rec.Seq); i += lineLength {
			if i+lineLength > len(rec.Seq) {
				_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(rec.Seq[i:]))
				common.ExitIfError(err)
			} else {
				_, err = fmt.Fprintf(file, "%s\n", dna.BasesToString(rec.Seq[i:i+lineLength]))
				common.ExitIfError(err)
			}
		}
	}
}

func WriteToGraphHandle(file io.Writer, gg *SimpleGraph, lineLength int) {
	var err error
	var i, j int
	for i = 0; i < len(gg.Nodes); i++ {
		_, err = fmt.Fprintf(file, "%s\n", ">"+gg.Nodes[i].Name+":"+fmt.Sprint(gg.Nodes[i].Id))
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
		_, err = fmt.Fprintf(file, "%s", gg.Nodes[i].Name+":"+fmt.Sprint(gg.Nodes[i].Id))
		near := gg.Nodes[i].Next
		for j = 0; j < len(near); j++ {
			_, err = fmt.Fprintf(file, "\t%v\t%s", near[j].Prob, near[j].Dest.Name+":"+fmt.Sprint(near[j].Dest.Id))
			common.ExitIfError(err)
		}
		_, err = fmt.Fprintf(file, "\n")
		common.ExitIfError(err)
	}
}

func PrintGraph(gg *SimpleGraph) {
	lineLength := 50
	var i, j int
	var name string
	for i = 0; i < len(gg.Nodes); i++ {
		name = ">" + gg.Nodes[i].Name + ":" + fmt.Sprint(gg.Nodes[i].Id)
		fmt.Printf("%s\n", name)
		for j = 0; j < len(gg.Nodes[i].Seq); j += lineLength {
			if j+lineLength > len(gg.Nodes[i].Seq) {
				fmt.Printf("%s\n", dna.BasesToString(gg.Nodes[i].Seq[j:]))

			} else {
				fmt.Printf("%s\n", dna.BasesToString(gg.Nodes[i].Seq[j:j+lineLength]))
			}
		}
	}
	for i = 0; i < len(gg.Nodes); i++ {
		fmt.Printf("%s", gg.Nodes[i].Name+":"+fmt.Sprint(gg.Nodes[i].Id))
		near := gg.Nodes[i].Next
		for j = 0; j < len(near); j++ {
			fmt.Printf("\t%v\t%s", near[j].Prob, near[j].Dest.Name+":"+fmt.Sprint(near[j].Dest.Id))
		}
		fmt.Print("\n")
	}
}

func Write(filename string, sg *SimpleGraph) {
	lineLength := 50
	file := fileio.EasyCreate(filename)
	defer file.Close()
	WriteToGraphHandle(file, sg, lineLength)
}
