package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"sync"
)

type Node struct {
	qSeq *[]dna.Base
}

type GenomeGraph struct {
	nodes []*Node
	edges map[Node][]*Node
	lock  sync.RWMutex
}

func NewGraph() *GenomeGraph {
	graph := new(GenomeGraph)
	graph.nodes = make([]*Node, 0)
	graph.edges = make(map[Node][]*Node, 0)
	return graph
}

func (g *GenomeGraph) AddNode(n *Node) {
	g.lock.Lock()
	g.nodes = append(g.nodes, n)
	g.lock.Unlock()
}

func (g *GenomeGraph) AddEdge(u, v *Node) {
	g.lock.Lock()
	if g.edges == nil {
		g.edges = make(map[Node][]*Node)
	}
	g.edges[*u] = append(g.edges[*u], v)
	g.lock.Unlock()
}

func (n *Node) String() string {
	return fmt.Sprintf("%s", dna.BasesToString(*n.qSeq))
}

func (g *GenomeGraph) String() {
	g.lock.RLock()
	s := ""
	for i := 0; i < len(g.nodes); i++ {
		s += g.nodes[i].String()
		near := g.edges[*g.nodes[i]]
		for j := 0; j < len(near); j++ {
			s += " -> " + near[j].String()
		}
		s += "\n"
	}
	fmt.Print(s)
	g.lock.RUnlock()
}

func FillGraph(input [][]*fasta.Fasta) *GenomeGraph {
	g := NewGraph()
	var prev *Node
	var curr *Node
	for i := 0; i < len(input); i++ {
		for j := 0; j < len(input[i]); j++ {
			curr = &Node{qSeq: &input[i][j].Seq}
			g.AddNode(curr)
			if j > 0 {
				g.AddEdge(prev, curr)
			}
			prev = curr
		}
	}
	return g
}

func usage() {
	fmt.Print(
		"GenomeGraph dev - b\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	inFile := flag.Arg(0)
	records, _ := fasta.ReadNew(inFile)
	a := fasta.DivideFastaAll(records, 3)
	g := FillGraph(a)
	g.String()
	//newGenomeGraph(inFile)
}
/*
func newGenomeGraph(inFile string) {
	records, err := fasta.ReadNew(inFile)

	for j, _ := range records {
		fmt.Println(dna.BasesToString(records[j].Seq))

	}
	a := fasta.DivideFastaAll(records, 3)
	g := FillGraph(a)
	for i, _ := range a {
		for j, _ := range a[i] {
			fmt.Println(dna.BasesToString(a[i][j].Seq))
		}
		fmt.Printf("\n")
	}
	if err != nil {
		log.Fatal(err)
	}
	g.String()
}
*/
