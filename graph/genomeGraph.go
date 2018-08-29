package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/qDna"
	"log"
	"sync"
)

type Node struct {
	qSeq *qDna.QFrag
}

type Edge struct {
	Next *Node
	Prob float64
}

type GenomeGraph struct {
	nodes   []*Node
	edges   map[*Node][]*Edge
	numNode int64
	lock    sync.RWMutex
}

func NewGraph() *GenomeGraph {
	graph := new(GenomeGraph)
	graph.nodes = make([]*Node, 0)
	graph.edges = make(map[*Node][]*Edge, 0)
	return graph
}

func (g *GenomeGraph) AddNode(n *Node) {
	g.lock.Lock()
	g.nodes = append(g.nodes, n)
	g.numNode = g.numNode + 1
	g.lock.Unlock()
}

func (g *GenomeGraph) AddEdge(u, v *Node, p float64) {
	g.lock.Lock()
	if g.edges == nil {
		g.edges = make(map[*Node][]*Edge)
	}
	tmp := Edge{Next: v, Prob: p}
	g.edges[u] = append(g.edges[u], &tmp)
	g.lock.Unlock()
}

func (n *Node) String() string {
	//return fmt.Sprintf("%v", dna.BasesToString(qDna.ToFasta(n.qSeq).Seq))
	var s string
	for i := 0; i < len(n.qSeq.Seq); i++ {
		s += fmt.Sprintf("%v", n.qSeq.Seq[i])
		s += ", "
	}
	return s
}

func (g *GenomeGraph) String() string {
	g.lock.RLock()
	s := ""
	for i := 0; i < len(g.nodes); i++ {
		s += g.nodes[i].String()
		near := g.edges[g.nodes[i]]
		for j := 0; j < len(near); j++ {
			s += " -> " + near[j].Next.String() + "\n"
		}
		s += "\n"
	}

	//fmt.Print(s)
	g.lock.RUnlock()
	return s
}

func FillGraph(input [][]*fasta.Fasta) *GenomeGraph {
	g := NewGraph()
	var prev *Node
	var curr *Node
	for i := 0; i < len(input); i++ {
		for j := 0; j < len(input[i]); j++ {
			curr = &Node{qSeq: qDna.FromFasta(input[i][j])}
			g.AddNode(curr)
			if j > 0 {
				g.AddEdge(prev, curr, 1)
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
	//records, _ := fasta.ReadNew(inFile)
	//a := fasta.DivideFastaAll(records, 3)
	//g := FillGraph(a)
	//g.String()
	newGenomeGraph(inFile)
}

//Modified from align package
func AlignToGraph(alpha *fasta.Fasta, beta *fasta.Fasta) *GenomeGraph {
	g := NewGraph()
	var curr *Node
	var prev *Node
	var lastMatch *Node = nil
	_, route := align.AffineGap(alpha.Seq, beta.Seq, align.DefaultScoreMatrix, -400, -30)
	var pairSeq []*fasta.Fasta
	pairSeq = append(pairSeq, alpha)
	pairSeq = append(pairSeq, beta)
	mAlign := align.AllSeqAffine(pairSeq, align.DefaultScoreMatrix, -400, -30)
	var seqIdx int64 = 0
	var prevOper align.ColType
	for i, _ := range route {
		oper := route[i].Op
		switch oper {
		case align.ColM:
			a := qDna.FromDnaToQFrag(mAlign[0].Seq[seqIdx:seqIdx+route[i].RunLength], mAlign[0].Name)
			b := qDna.FromDnaToQFrag(mAlign[1].Seq[seqIdx:seqIdx+route[i].RunLength], mAlign[1].Name)
			curr = &Node{qSeq: qDna.PairwiseAverage(a, b)}
			g.AddNode(curr)
			if g.numNode > 0 {
				g.AddEdge(prev, curr, 1)
			}
			if prevOper == align.ColI || prevOper == align.ColD {
				g.AddEdge(lastMatch, curr, 1)
			}
			//For debugging purposes:
			//tmp := qDna.PairwiseAverage(a, b)
			//for i, _ := range tmp.Seq {
			//	fmt.Println(tmp.Seq[i])
			//}
			//fmt.Println("end")
			lastMatch = curr
		//case 1:
		case align.ColI:
			curr = &Node{qSeq: qDna.FromDnaToQFrag(mAlign[1].Seq[seqIdx:seqIdx+route[i].RunLength], mAlign[1].Name)}
			g.AddNode(curr)
			if g.numNode > 0 {
				g.AddEdge(prev, curr, 1)
			}
			//case 2:
		case align.ColD:
			curr = &Node{qSeq: qDna.FromDnaToQFrag(mAlign[0].Seq[seqIdx:seqIdx+route[i].RunLength], mAlign[0].Name)}
			g.AddNode(curr)
			if g.numNode > 0 {
				g.AddEdge(lastMatch, curr, 1)
			}
		}
		seqIdx = seqIdx + route[i].RunLength
		prev = curr
		prevOper = oper
	}
	//fmt.Println(route)
	fmt.Println("Input Sequences Aligned:")
	fmt.Println(dna.BasesToString(mAlign[0].Seq))
	fmt.Println(dna.BasesToString(mAlign[1].Seq))
	return g
}

func newGenomeGraph(inFile string) {
	records, err := fasta.ReadNew(inFile)
	fmt.Println("Input Sequences: ")
	for j, _ := range records {
		fmt.Println(dna.BasesToString(records[j].Seq))
	}
	//a := fasta.DivideFastaAll(records, 3)
	//g := FillGraph(a)
	if err != nil {
		log.Fatal(err)
	}
	//s := g.String()
	//fmt.Println(s)
	//mAlign := align.AllSeqAffine(records, align.DefaultScoreMatrix, -400, -30)
	alpha := records[0]
	beta := records[1]
	//z, cigar := align.AffineGap(alpha, beta, align.DefaultScoreMatrix, -400, -30)
	//fmt.Println(dna.BasesToString(mAlign[0].Seq))
	//fmt.Println(dna.BasesToString(mAlign[1].Seq))
	//fmt.Println(z)
	//fmt.Println(cigar)
	gg := AlignToGraph(alpha, beta)
	fmt.Println(gg.String())

}
