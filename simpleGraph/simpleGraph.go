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
	Nodes   []*Node
	//Edges   map[*Node][]*Edge
	//lock    sync.RWMutex
}

type Node struct {
	Id  int64
	Name string
	Seq []dna.Base
	Prev []*Edge
	Next []*Edge
}

type Edge struct {
	//Curr *Node
	Next *Node
	Prob float64
	//Strand bool
}

func AddNode(g *SimpleGraph, n *Node) {
	//g.lock.Lock()
	g.Nodes = append(g.Nodes, n)
	//g.lock.Unlock()
}
//func AddEdge(u, v *Node, uStrand, vStrand bool, p float64) {
func AddEdge(u, v *Node, p float64) {
	//g.lock.Lock()
	//if g.Edges == nil {
	//	g.Edges = make(map[*Node][]*Edge)
	//}
	//if uStand == true {
	//	u.Next = append(u.Next, &Edge{Next: v, Prob: p, Strand: vStand})
	//}
	//if uStand == false {
	//	u.Next = append(u.Next, &Edge{Next: v, Prob: p, Strand: uStand})
	//}
	//if vStand == true {

	//}
	//v.Prev = append(g.Edges[v], &Edge{Next: u, Prob: p, Strand: false})
	//g.lock.Unlock()
	u.Next = append(u.Next, &Edge{Next: v, Prob: p})
	v.Prev = append(v.Prev, &Edge{Next: u, Prob: p})
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
			if seqIdx > 0 {
				AddEdge(answer.Nodes[seqIdx-1], &tmp, 1)
			}
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
	var startBase int64 = 1
	lineLength := 50
	var i, j int
	for i = 0; i < len(gg.Nodes); i++ {
		fmt.Printf("%s%s\t%v\t%v\n", ">", gg.Nodes[i].Name, startBase, startBase+int64(len(gg.Nodes[i].Seq)))
		startBase += int64(len(gg.Nodes[i].Seq)) + 1
		for j = 0; j < len(gg.Nodes[i].Seq);j+=lineLength {
			if j+lineLength > len(gg.Nodes[i].Seq) {
				fmt.Printf("%s\n", dna.BasesToString(gg.Nodes[i].Seq[j:]))
			} else {
				fmt.Printf("%s\n", dna.BasesToString(gg.Nodes[i].Seq[j:j+lineLength]))
			}
		}
	}
	for i = 0; i < len(gg.Nodes); i++ {
		fmt.Printf("%s\t", gg.Nodes[i].Name)
		//near := gg.Edges[gg.Nodes[i]]
		near := gg.Nodes[i].Next
		for j = 0; j < len(near);j++ {
			
			fmt.Printf("%s\t", near[j].Next.Name)
			
			
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
