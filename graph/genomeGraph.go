package graph

import (
	"bufio"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/qDna"
	"github.com/vertgenlab/gonomics/vcf"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
	"sync"
)

type Node struct {
	QSeq *qDna.QFrag
	//unique identifier
	NodeID int64
}

type Edge struct {
	Next *Node
	Prob float64
}

type GenomeGraph struct {
	Nodes   []*Node
	Edges   map[*Node][]*Edge
	NumNode int64
	lock    sync.RWMutex
}

func NewGraph() *GenomeGraph {
	graph := new(GenomeGraph)
	graph.Nodes = make([]*Node, 0)
	graph.Edges = make(map[*Node][]*Edge, 0)
	return graph
}

func (g *GenomeGraph) AddNode(n *Node) {
	g.lock.Lock()
	g.Nodes = append(g.Nodes, n)
	g.NumNode = g.NumNode + 1
	g.lock.Unlock()
}

func (g *GenomeGraph) AddEdge(u, v *Node, p float64) {
	g.lock.Lock()
	if g.Edges == nil {
		g.Edges = make(map[*Node][]*Edge)
	}
	tmp := Edge{Next: v, Prob: p}
	g.Edges[u] = append(g.Edges[u], &tmp)
	g.lock.Unlock()
}

func SnpToSeq(sequence []dna.Base, snps []*vcf.Vcf) []dna.Base {
	answer := []dna.Base{}
	for s := 0; s < len(sequence); s++ {
		answer = append(answer, sequence[s])
	}
	for i := 0; i < len(snps); i++ {
		base := dna.StringToBases(snps[i].Alt[0])
		answer[snps[i].Pos-1] = base[0]
	}
	return answer
}

func SnpQFrag(alpha []dna.Base, beta []dna.Base, name string, start int64, end int64) *qDna.QFrag {
	return qDna.PairwiseAverage(qDna.FromDnaToQFrag(alpha, name), qDna.FromDnaToQFrag(beta, name), start, end, name)
}

func SeqToGraph(vcfFile []*vcf.Vcf, sequence *fasta.Fasta, gsw *GenomeGraph) *GenomeGraph {
	g := gsw
	var curr *Node
	var currMatch *Node
	var prev *Node
	var lastMatch *Node = nil
	var lastPos int64 = 0
	var snps []*vcf.Vcf
	lastPos = 0
	for i := 0; i < len(vcfFile); i++ {
		if strings.Contains(vcfFile[i].Info, "SVTYPE=SNP") {
			snps = append(snps, vcfFile[i])
		}
		//case insertion
		if strings.Contains(vcfFile[i].Info, "SVTYPE=INS") {
			//logic for calling SNPs
			snpSeq := SnpToSeq(sequence.Seq, snps)
			//added location of fragment and converted to 1 base by add 1
			currMatch = &Node{QSeq: SnpQFrag(sequence.Seq[lastPos:vcfFile[i].Pos], snpSeq[lastPos:vcfFile[i].Pos], sequence.Name, lastPos+1, int64(vcfFile[i].Pos)), NodeID: g.NumNode}
			snps = nil
			g.AddNode(currMatch)
			bases := dna.StringToBases(vcfFile[i].Alt[0])
			curr = &Node{QSeq: qDna.QFragCoord(bases[1:], sequence.Name, lastPos+1, lastPos+1), NodeID: g.NumNode}
			g.AddNode(curr)
			lastPos++
			if lastMatch != nil {
				g.AddEdge(lastMatch, currMatch, 1)
				g.AddEdge(prev, currMatch, 1)
			}
			g.AddEdge(currMatch, curr, 1)
			lastMatch = currMatch
		}
		//deletion in vcf record
		if strings.Contains(vcfFile[i].Info, "SVTYPE=DEL") {

			snpSeq := SnpToSeq(sequence.Seq, snps)
			currMatch = &Node{QSeq: SnpQFrag(sequence.Seq[lastPos:vcfFile[i].Pos], snpSeq[lastPos:vcfFile[i].Pos], sequence.Name, lastPos+1, int64(vcfFile[i].Pos)), NodeID: g.NumNode}
			snps = nil
			g.AddNode(currMatch)

			bases := dna.StringToBases(vcfFile[i].Ref)
			curr = &Node{QSeq: qDna.QFragCoord(bases[1:], sequence.Name, int64(vcfFile[i].Pos+1), int64(vcfFile[i].Pos)+int64(len(bases))-1), NodeID: g.NumNode}
			g.AddNode(curr)
			lastPos = int64(vcfFile[i].Pos) + int64(len(bases)) - 1
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
	lastNode := &Node{QSeq: SnpQFrag(sequence.Seq[lastPos:len(sequence.Seq)], snpSeq[lastPos:len(sequence.Seq)], sequence.Name, lastPos+1, int64(len(sequence.Seq))), NodeID: g.NumNode}
	g.AddNode(lastNode)
	g.AddEdge(prev, lastNode, 1)
	g.AddEdge(lastMatch, lastNode, 1)
	return g
}

func RefernceToGraph(vcfFile []*vcf.Vcf, reference []*fasta.Fasta, g *GenomeGraph) *GenomeGraph {
	gsw := g
	vcfSplit := vcf.VcfSplit(vcfFile, reference)
	if len(vcfSplit) != len(reference) {
		fmt.Println("Slice of vcfs do not equal reference length")
	}
	for i := 0; i < len(reference); i++ {
		gsw = SeqToGraph(vcfSplit[i], reference[i], gsw)
	}
	return gsw
}

//Dictionary look up of sequence by name
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
	//Use Bases to string to visualize qDNA, does not work with qDNA are SNPs
	return fmt.Sprintf("%v", dna.BasesToString(qDna.ToFasta(n.QSeq).Seq))
	//var s string
	//for i := 0; i < len(n.QSeq.Seq); i++ {
	//	s += fmt.Sprintf("%v", n.QSeq.Seq[i])
	//	s += ", "
	//}
	//return s
}

func (g *GenomeGraph) String() string {
	g.lock.RLock()
	s := ""
	for i := 0; i < len(g.Nodes); i++ {
		s += g.Nodes[i].String()
		near := g.Edges[g.Nodes[i]]
		for j := 0; j < len(near); j++ {
			s += " -> " + near[j].Next.String()
		}
		s += "\n"
	}
	g.lock.RUnlock()
	return s
}

func PrintGraph(g *GenomeGraph) {
	for i := 0; i < len(g.Nodes); i++ {
		id := strconv.FormatInt(g.Nodes[i].NodeID, 10)
		fmt.Printf("%s %s %v %v\n", "> "+id, g.Nodes[i].QSeq.From[0].Chr, g.Nodes[i].QSeq.From[0].Start, g.Nodes[i].QSeq.From[0].End)
		for j := 0; j < len(g.Nodes[i].QSeq.Seq); j++ {
			fmt.Printf("%v %v %v %v\n", g.Nodes[i].QSeq.Seq[j].A, g.Nodes[i].QSeq.Seq[j].C, g.Nodes[i].QSeq.Seq[j].G, g.Nodes[i].QSeq.Seq[j].T)
		}
		fmt.Print("\n")
	}
	for x := 0; x < len(g.Nodes); x++ {
		fmt.Printf("%v", g.Nodes[x].NodeID)
		near := g.Edges[g.Nodes[x]]
		for y := 0; y < len(near); y++ {
			fmt.Printf(" %v", near[y].Next.NodeID)
		}
		fmt.Print("\n")
	}
}

func PrintGraphRecord(gsw []*GenomeGraph) {
	for i := 0; i < len(gsw); i++ {
		PrintGraph(gsw[i])
	}
}

func WriteGraphToFileHandle(file *os.File, input *GenomeGraph) error {
	var err error

	for i := 0; i < len(input.Nodes); i++ {
		id := strconv.FormatInt(input.Nodes[i].NodeID, 10)
		_, err = fmt.Fprintf(file, "%s %s %v %v\n", "> "+id, input.Nodes[i].QSeq.From[0].Chr, input.Nodes[i].QSeq.From[0].Start, input.Nodes[i].QSeq.From[0].End)
		for j := 0; j < len(input.Nodes[i].QSeq.Seq); j++ {
			_, err = fmt.Fprintf(file, "%v %v %v %v\n", input.Nodes[i].QSeq.Seq[j].A, input.Nodes[i].QSeq.Seq[j].C, input.Nodes[i].QSeq.Seq[j].G, input.Nodes[i].QSeq.Seq[j].T)
		}
		_, err = fmt.Fprintf(file, "\n")
	}
	for x := 0; x < len(input.Nodes); x++ {
		_, err = fmt.Fprintf(file, "%v", input.Nodes[x].NodeID)
		near := input.Edges[input.Nodes[x]]
		for y := 0; y < len(near); y++ {
			_, err = fmt.Fprintf(file, " %v", near[y].Next.NodeID)
		}
		_, err = fmt.Fprintf(file, "\n")
	}
	common.ExitIfError(err)
	return err
}

func Write(filename string, graph *GenomeGraph) {
	file, err := os.Create(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	WriteGraphToFileHandle(file, graph)
}

func Read(filename string) *GenomeGraph {
	var answer *GenomeGraph
	m := make(map[int64]*Node)
	file, _ := os.Open(filename)
	defer file.Close()
	reader := bufio.NewReader(file)
	answer = NewGraph()
	var err error
	var line string
	var currNode *Node
	var words []byte
	for ; err != io.EOF; words, _, err = reader.ReadLine() {
		line = string(words[:])
		data := strings.Split(line, " ")
		if strings.HasPrefix(line, ">") {
			//New Node
			data = strings.Split(line, " ")
			start, _ := strconv.ParseInt(data[3], 10, 64)
			end, _ := strconv.ParseInt(data[4], 10, 64)
			loc := qDna.Location{Assembly: "", Chr: data[2], Start: start, End: end}
			currNode = &Node{QSeq: &qDna.QFrag{Seq: nil, From: []*qDna.Location{&loc}, Fwd: nil, Rev: nil}, NodeID: answer.NumNode}
			answer.AddNode(currNode)
			_, ok := m[currNode.NodeID]
			if !ok {
				m[currNode.NodeID] = currNode
			}
		} else if len(data) == 4 {
			data = strings.Split(line, " ")
			a, _ := strconv.ParseFloat(data[0], 64)
			c, _ := strconv.ParseFloat(data[1], 64)
			g, _ := strconv.ParseFloat(data[2], 64)
			t, _ := strconv.ParseFloat(data[3], 64)
			currNode.QSeq.Seq = append(currNode.QSeq.Seq, &qDna.QBase{A: a, C: c, G: g, T: t})
		} else if len(data) == 3 {
			e1, _ := strconv.ParseInt(data[0], 10, 64)
			e2, _ := strconv.ParseInt(data[1], 10, 64)
			e3, _ := strconv.ParseInt(data[2], 10, 64)
			answer.AddEdge(m[e1], m[e2], 1)
			answer.AddEdge(m[e1], m[e3], 1)

		} else if len(data) == 2 {
			e1, _ := strconv.ParseInt(data[0], 10, 64)
			e2, _ := strconv.ParseInt(data[1], 10, 64)
			answer.AddEdge(m[e1], m[e2], 1)
		} else {

		}
	}
	return answer
}
