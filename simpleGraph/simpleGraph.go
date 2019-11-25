package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"strings"
	//"sync"
)

type SimpleGraph struct {
	Nodes []*Node
	//Edges   map[*Node][]*Edge
	//lock    sync.RWMutex
}

type Node struct {
	Id   int64
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
	//g.lock.Lock()
	g.Nodes = append(g.Nodes, n)
	//g.lock.Unlock()
}

//func AddEdge(u, v *Node, uStrand, vStrand bool, p float64) {
func AddEdge(u, v *Node, p float64) {
	//g.lock.Lock()
	//g.lock.Unlock()
	u.Next = append(u.Next, &Edge{Dest: v, Prob: p})
	v.Prev = append(v.Prev, &Edge{Dest: u, Prob: p})
}

//func chaining()

//func Read(filename string) ([]*Node, []string) {
func Read(filename string) *SimpleGraph {
	answer := NewGraph()
	var line string
	var currSeq []dna.Base
	//var answer []*Node
	var seqIdx int64 = -1
	var doneReading bool = false
	//var chromIdx int64 = 0
	file := fileio.EasyOpen(filename)
	defer file.Close()
	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, ">") {
			seqIdx++
			//tmp := Node{Id: common.StringToInt64(line[1:len(line)])}
			tmp := Node{Id: seqIdx, Name: line[1:len(line)]}
			//answer = append(answer, &tmp)
			AddNode(answer, &tmp)
			//if seqIdx > 0 {
			//	AddEdge(answer.Nodes[seqIdx-1], &tmp, 1)
			//}
		} else {
			currSeq = dna.StringToBases(line)
			dna.AllToUpper(currSeq)
			answer.Nodes[seqIdx].Seq = append(answer.Nodes[seqIdx].Seq, currSeq...)
		}
	}
	return answer
}

func NewGraph() *SimpleGraph {
	graph := new(SimpleGraph)
	graph.Nodes = make([]*Node, 0)
	//graph.Edges = make(map[*Node][]*Edge, 0)
	return graph
}

/*
func AddEdgeTmp(u, v *Node) {
	u.Next, v.Prev = append(u.Next, v), append(v.Prev, u)
}*/

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

func PrintGraph(gg *SimpleGraph) {

	lineLength := 50
	var i, j int
	for i = 0; i < len(gg.Nodes); i++ {
		fmt.Printf("%s%s\n", ">", gg.Nodes[i].Name)
		for j = 0; j < len(gg.Nodes[i].Seq); j += lineLength {
			if j+lineLength > len(gg.Nodes[i].Seq) {
				fmt.Printf("%s\n", dna.BasesToString(gg.Nodes[i].Seq[j:]))
			} else {
				fmt.Printf("%s\n", dna.BasesToString(gg.Nodes[i].Seq[j:j+lineLength]))
			}
		}
	}
	for i = 0; i < len(gg.Nodes); i++ {
		fmt.Printf("%s\t", gg.Nodes[i].Name)
		near := gg.Nodes[i].Next
		for j = 0; j < len(near); j++ {
			fmt.Printf("%s\t", near[j].Dest.Name)
		}
		fmt.Print("\n")
	}
}

func Write(filename string, records []*Node) {
	lineLength := 50
	file := fileio.EasyCreate(filename)
	defer file.Close()

	WriteToFileHandle(file, records, lineLength)
}
