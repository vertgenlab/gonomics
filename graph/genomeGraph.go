package graph

import (
	"fmt"
	"sync"
	"strings"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/qDna"
	"github.com/vertgenlab/gonomics/vcf"
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

func SnpToSeq(sequence []dna.Base, snps []*vcf.Vcf) []dna.Base{
	answer := []dna.Base{}
	for s := 0; s < len(sequence); s++ {
		answer = append(answer,sequence[s])
	}
	for i := 0; i < len(snps); i++ {
		base := dna.StringToBases(snps[i].Alt)
		answer[snps[i].Pos-1] = base[0]
	}
	return answer
}

func SnpQFrag(alpha []dna.Base, beta []dna.Base, name string) *qDna.QFrag{
	return qDna.PairwiseAverage(qDna.FromDnaToQFrag(alpha, name), qDna.FromDnaToQFrag(beta, name))
}

func SeqToGraph(vcfFile []*vcf.Vcf, sequence *fasta.Fasta) *GenomeGraph {
	g := NewGraph()
	var curr *Node
	var currMatch *Node
	var prev *Node
	var lastMatch *Node = nil
	var lastPos int64
	var snps []*vcf.Vcf
	lastPos = 0
	for i := 0; i < len(vcfFile); i++ {
		if strings.Compare(vcfFile[i].Format, "SVTYPE=SNP") == 0 {
			snps = append(snps, vcfFile[i])
		}
		//case insertion
		if strings.Compare(vcfFile[i].Format, "SVTYPE=INS") == 0 {
			//logic for calling SNPs
			snpSeq := SnpToSeq(sequence.Seq, snps)
			currMatch = &Node{qSeq: SnpQFrag(sequence.Seq[lastPos:vcfFile[i].Pos], snpSeq[lastPos:vcfFile[i].Pos], sequence.Name)}
			snps = nil
			g.AddNode(currMatch)

			bases := dna.StringToBases(vcfFile[i].Alt)
			curr = &Node{qSeq: qDna.FromDnaToQFrag(bases[1:len(bases)], sequence.Name)}
			g.AddNode(curr)
			if lastMatch != nil {
				g.AddEdge(lastMatch, currMatch, 1)
				g.AddEdge(prev, currMatch, 1)
			}
			g.AddEdge(currMatch, curr, 1)
			lastMatch = currMatch
		}
		//deletion in vcf record
		if strings.Compare(vcfFile[i].Format, "SVTYPE=DEL") == 0 {
			
			snpSeq := SnpToSeq(sequence.Seq, snps)
			currMatch = &Node{qSeq: SnpQFrag(sequence.Seq[lastPos:vcfFile[i].Pos], snpSeq[lastPos:vcfFile[i].Pos], sequence.Name)}
			snps = nil
			g.AddNode(currMatch)

			bases := dna.StringToBases(vcfFile[i].Ref)
			curr = &Node{qSeq: qDna.FromDnaToQFrag(bases[1:len(bases)], sequence.Name)}
			g.AddNode(curr)
			lastPos = vcfFile[i].Pos + int64(len(bases)) - 1
			if lastMatch != nil {
				g.AddEdge(lastMatch, currMatch, 1)
				g.AddEdge(prev, currMatch, 1)
			}
			g.AddEdge(currMatch, curr, 1)
			lastMatch = currMatch
		}
		prev = curr
	}

	//last match case:
	snpSeq := SnpToSeq(sequence.Seq, snps)
	lastNode := &Node{qSeq: SnpQFrag(sequence.Seq[lastPos:len(sequence.Seq)], snpSeq[lastPos:len(sequence.Seq)], sequence.Name)}
	g.AddNode(lastNode)
	g.AddEdge(prev, lastNode, 1)
	g.AddEdge(lastMatch, lastNode, 1)
	return g
}

func RefernceToGraph(vcfFile []*vcf.Vcf, reference []*fasta.Fasta) []*GenomeGraph {
	var gGraph []*GenomeGraph
	var curr *GenomeGraph
	vcfSplit := vcf.VcfSplit(vcfFile, reference)
	if len(vcfSplit) != len(reference) {
		fmt.Println("Slice of vcfs do not equal reference length")
	}
	for i := 0; i < len(reference); i++ {
		curr = SeqToGraph(vcfSplit[i], reference[i])
		gGraph = append(gGraph, curr)
	}
	return gGraph
}

func FastaMap(ref []*fasta.Fasta) map[string][]dna.Base {
	m := make(map[string][]dna.Base)
	//var answer []*fasta.Fasta
	var curr *fasta.Fasta
	for i := 0; i < len(ref); i++ {
		curr = ref[i]
		_, ok := m[curr.Name]
		if !ok {
			m[curr.Name] = curr.Seq
		}
	}
	return m
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
			s += " -> " + near[j].Next.String()
		}
		s += "\n"
	}
	g.lock.RUnlock()
	return s
}