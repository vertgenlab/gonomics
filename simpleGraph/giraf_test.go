package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/giraf"
	"log"
	"testing"
)

// TestGraph Structure
//             n2          e0 = 1
//         e1/    \e2      e1 = 0.05
//      e0  /  e3  \       e2 = 1
//  n0 --- n1 ----- n4     e3 = 0.8
//          \      /       e4 = 0.15
//         e4\    /e5      e5 = 1
//             n3
//
//               A
//             /    \
//            /      \
//  ATG --- CG ----- TAA
//            \      /
//             \    /
//               T

// Sam Header:
//@HD     VN:1.6  SO:coordinate
//@SQ     SN:n0 LN:3
//@SQ     SN:n1 LN:2
//@SQ     SN:n2 LN:1
//@SQ     SN:n3 LN:1
//@SQ     SN:n4 LN:3

// Test Functions
func MakeTestGraph() *SimpleGraph {
	graph := NewGraph()

	var n0, n1, n2, n3, n4 *Node
	var e0, e1, e2, e3, e4, e5 *Edge

	// Make Nodes
	n0 = &Node{
		Id:   0,
		Name: "n0",
		Seq:  dna.StringToBases("ATG")}

	n1 = &Node{
		Id:   1,
		Name: "n1",
		Seq:  dna.StringToBases("CG")}

	n2 = &Node{
		Id:   2,
		Name: "n2",
		Seq:  dna.StringToBases("A")}

	n3 = &Node{
		Id:   3,
		Name: "n3",
		Seq:  dna.StringToBases("T")}

	n4 = &Node{
		Id:   4,
		Name: "n4",
		Seq:  dna.StringToBases("TAA")}

	// Make Edges
	e0 = &Edge{
		Dest: n1,
		Prob: 1}

	e1 = &Edge{
		Dest: n2,
		Prob: 0.05}

	e2 = &Edge{
		Dest: n4,
		Prob: 1}

	e3 = &Edge{
		Dest: n4,
		Prob: 0.8}

	e4 = &Edge{
		Dest: n3,
		Prob: 0.15}

	e5 = &Edge{
		Dest: n4,
		Prob: 1}

	// Define Paths
	n0.Next = append(n0.Next, e0)
	n1.Next = append(n1.Next, e1, e3, e4)
	n1.Prev = append(n1.Prev, e0)
	n2.Next = append(n2.Next, e2)
	n2.Prev = append(n2.Prev, e1)
	n3.Next = append(n3.Next, e5)
	n3.Prev = append(n3.Prev, e4)
	n4.Prev = append(n4.Prev, e2, e3, e5)

	graph.Nodes = append(graph.Nodes, n0, n1, n2, n3, n4)

	for i := 0; i < len(graph.Nodes); i++ {
		graph.Nodes[i].SeqTwoBit = dnaTwoBit.NewTwoBit(graph.Nodes[i].Seq)
	}

	return graph
}

// check struct generated with parameters reads := RandGiraf(MakeTestGraph(), 1, 4, seed)
var check = giraf.Giraf{
	QName:     "0_3_2_2_-",
	QStart:    0,
	QEnd:      4,
	PosStrand: false,                                                     // rev strand, must reverse complement
	Path:      &giraf.Path{TStart: 2, Nodes: []uint32{0, 1, 2}, TEnd: 1}, // Nodes 0->1->2, start base 3, end base 1
	Aln:       []*cigar.Cigar{{RunLength: 4, Op: 'M', Sequence: nil}},
	AlnScore:  16607,
	MapQ:      30,
	Seq:       []dna.Base{3, 1, 2, 1}, // TCGC
	Qual:      []uint8{16, 38, 38, 36},
	Notes:     nil}

func TestRandGiraf(t *testing.T) {
	var seed int64 = 777

	// Uncomment following line for random reads
	//seed = time.Now().UnixNano()

	reads := RandGiraf(MakeTestGraph(), 10, 4, seed)

	if reads[0].QName != check.QName {
		log.Fatalln("Reads do not match")
	}
	if reads[0].PosStrand != check.PosStrand {
		log.Fatalln("Reads do not match")
	}
	if reads[0].Path.TStart != check.Path.TStart ||
		reads[0].Path.TEnd != check.Path.TEnd ||
		reads[0].Path.Nodes[0] != check.Path.Nodes[0] ||
		reads[0].Path.Nodes[len(reads[0].Path.Nodes)-1] != check.Path.Nodes[len(reads[0].Path.Nodes)-1] {
		log.Fatalln("Reads do not match")
	}
	if dna.CompareSeqsIgnoreCase(reads[0].Seq, check.Seq) == 1 {
		log.Fatalln("Reads do not match")
	}

	fmt.Println("Sequences Match")

	RandSomaticMutations(MakeTestGraph(), reads, 1, 0.75, seed)

	if reads[0].Seq[1] != 1 {
		log.Fatalln("Problem with simulating somatic mutations")
	}

	fmt.Println("Somatic Mutations Generated Correctly")

	/*
		for i := 0; i < len(reads); i++ {
			fmt.Println(reads[i])
			fmt.Println(reads[i].Path, reads[i].Aln[0])
		}
	*/

}
