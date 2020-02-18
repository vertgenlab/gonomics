package alleles

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"math/rand"
	"testing"
	"time"
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

// Test Functions
func MakeTestGraph() *simpleGraph.SimpleGraph {
	graph := simpleGraph.NewGraph()

	var n0, n1, n2, n3, n4 *simpleGraph.Node
	var e0, e1, e2, e3, e4, e5 *simpleGraph.Edge

	// Make Nodes
	n0 = &simpleGraph.Node{
		Id: 	0,
		Name: 	"n0",
		Seq:	dna.StringToBases("ATG")}

	n1 = &simpleGraph.Node{
		Id: 	1,
		Name: 	"n1",
		Seq:	dna.StringToBases("CG")}

	n2 = &simpleGraph.Node{
		Id: 	2,
		Name: 	"n2",
		Seq:	dna.StringToBases("A")}

	n3 = &simpleGraph.Node{
		Id: 	3,
		Name: 	"n3",
		Seq:	dna.StringToBases("T")}

	n4 = &simpleGraph.Node{
		Id: 	4,
		Name: 	"n4",
		Seq:	dna.StringToBases("TAA")}

	// Make Edges
	e0 = &simpleGraph.Edge{
		Dest: 	n1,
		Prob: 	1}

	e1 = &simpleGraph.Edge{
		Dest: 	n2,
		Prob: 	0.05}

	e2 = &simpleGraph.Edge{
		Dest: 	n4,
		Prob: 	1}

	e3 = &simpleGraph.Edge{
		Dest: 	n4,
		Prob: 	0.8}

	e4 = &simpleGraph.Edge{
		Dest: 	n3,
		Prob: 	0.15}

	e5 = &simpleGraph.Edge{
		Dest: 	n4,
		Prob: 	1}

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

	return graph
}

func MakeTestAln() []*sam.SamAln {
	var reads = make([]*sam.SamAln, 0)
	var r1, r2, r3, r4, r5, r6, r7, r8, r9, r10 *sam.SamAln

	r1 = &sam.SamAln{
		RName: 	"n0",
		Seq:	dna.StringToBases("ATGCGTAA"),
		Extra: 	"0:1:4"}

	r2 = &sam.SamAln{
		RName: 	"n1",
		Seq:	dna.StringToBases("CGTAA"),
		Extra: 	"1:4"}

	r3 = &sam.SamAln{
		RName: 	"n1",
		Seq:	dna.StringToBases("CGT"),
		Extra: 	"1:3"}

	r4 = &sam.SamAln{
		RName: 	"n0",
		Seq:	dna.StringToBases("ATGCGA"),
		Extra: 	"0:1:2"}

	r5 = &sam.SamAln{
		RName: 	"n0",
		Seq:	dna.StringToBases("ATGCGTAA"),
		Extra: 	"0:1:4"}

	// Mutation in n1
	r6 = &sam.SamAln{
		RName: 	"n0",
		Seq:	dna.StringToBases("AAGCGTAA"),
		Extra: 	"0:1:4"}

	r7 = &sam.SamAln{
		RName: 	"n4",
		Seq:	dna.StringToBases("TAA"),
		Extra: 	"4"}

	r8 = &sam.SamAln{
		RName: 	"n3",
		Seq:	dna.StringToBases("A"),
		Extra: 	"3"}

	// Mutation in n5
	r9 = &sam.SamAln{
		RName: 	"n1",
		Seq:	dna.StringToBases("CGTTAG"),
		Extra: 	"1:3:4"}

	r10 = &sam.SamAln{
		RName: 	"n3",
		Seq:	dna.StringToBases("TTAA"),
		Extra: 	"3:4"}

	reads = append(reads, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
	return reads
}

func MakeTest() (*simpleGraph.SimpleGraph, []*sam.SamAln) {
	return MakeTestGraph(), MakeTestAln()
}

// Get random path through graph
func RandPath(graph *simpleGraph.SimpleGraph) []*simpleGraph.Node {
	answer := make([]*simpleGraph.Node, 0)
	rand.Seed(time.Now().UnixNano())
	var currentNode, nextNode *simpleGraph.Node

	currentNode = graph.Nodes[0]
	for {
		answer = append(answer, currentNode)

		if len(currentNode.Next) == 0 {
			break
		}
		nextNode = currentNode.Next[rand.Intn(len(currentNode.Next))].Dest
		currentNode = nextNode
	}
	return answer
}

// Currently does not make indels that are not in graph
// TODO: add in functionality for somatic SNVs
type NewSNV struct {
	Node 	*simpleGraph.Node
	Pos 	int
	Base 	dna.Base
}

func randBase (base dna.Base) dna.Base {
	rand.Seed(time.Now().UnixNano())

	var randomBase dna.Base

	for {
		randNum := rand.Intn(4)
		switch randNum {
		case 0 :
			randomBase = dna.A
		case 1 :
			randomBase = dna.C
		case 2 :
			randomBase = dna.G
		case 3 :
			randomBase = dna.T
		}

		if base != randomBase {
			break
		}
	}

	return randomBase
}

func MakeTestReads(graph *simpleGraph.SimpleGraph, numReads int, readLen int, numSNP int, numNewSNV int) []*sam.SamAln {
	answer := make([]*sam.SamAln, 0)
	rand.Seed(time.Now().UnixNano())

	var i, j int
	var possibleSNPids = make([]int, 0)
	var SNPids = make([]int, 0)
	// Choose random path through graph
	basePath := RandPath(graph)



	var NewSNVs []NewSNV

	for i = 0; i < numNewSNV; i++ {
		randNode := basePath[rand.Intn(len(basePath))]
		randPos := rand.Intn(len(randNode.Seq))
		randSub := randBase(randNode.Seq[randPos])

		NewSNVs = append(NewSNVs, NewSNV{randNode, randPos, randSub})
	}

	// Determine which nodes could diverge into a SNP
	for i = 0; i < len(basePath); i++ {
		if len(basePath[i].Next) > 1 {
			possibleSNPids = append(possibleSNPids, i)
		}
	}

	// Choose SNPs
	for i = 0; i < numSNP; i++ {
		randid := rand.Intn(len(possibleSNPids))
		var match bool = false
		for j = 0; j < len(SNPids); j++ {
			if possibleSNPids[randid] == SNPids[j] {
				match = true
				break
			}
		}
		if match == false {
			SNPids = append(SNPids, possibleSNPids[randid])
		}
	}

	// Pick a start node, if it is a snp then choose between two options
	var current *sam.SamAln
	for i = 0; i < numReads; i++ {
		startNodeId := rand.Intn(len(basePath))
		var altnode *simpleGraph.Node = nil

		for j = 0; j < len(SNPids); j++ {
			if startNodeId - 1 == SNPids[j] {
				for {
					randid := rand.Intn(len(basePath[SNPids[j]].Next))
					testSNPnode := basePath[SNPids[j]].Next[randid].Dest
					if testSNPnode != basePath[startNodeId] {
						altnode = testSNPnode
						break
					}
				}
				break
			}
		}
		randStart := rand.Intn(2)

		var readPath = make([]*simpleGraph.Node, 0)

		if altnode != nil && randStart == 1 {
			readPath = append(readPath, altnode)
		} else {
			readPath = append(readPath, basePath[startNodeId])
		}
		readStartPos := rand.Intn(len(readPath[0].Seq))

		var readSeq = make([]dna.Base, 0)
		readSeq = append(readSeq, readPath[0].Seq[readStartPos:]...)

		var currNodeId = startNodeId
		for len(readSeq) < readLen {

			if len(basePath[currNodeId].Next) == 0 {
				break
			}
			// Figure out if current node branches into a SNP
			altnode = nil
			for j = 0; j < len(SNPids); j++ {
				if currNodeId == SNPids[j] {
					for {
						testSNPnode := basePath[SNPids[j]].Next[rand.Intn(len(basePath[SNPids[j]].Next))].Dest
						if testSNPnode != basePath[currNodeId] {
							altnode = testSNPnode
							break
						}
					}
					break
				}
			}

			currNodeId++

			randStart = rand.Intn(2)
			if altnode != nil && randStart == 1 {
				readPath = append(readPath, altnode)
			} else {
				readPath = append(readPath, basePath[currNodeId])
			}


			missingSeq := readLen - len(readSeq)

			if len(basePath[currNodeId].Seq) > missingSeq {
				readSeq = append(readSeq, basePath[currNodeId].Seq[0:missingSeq]...)
			} else {
				readSeq = append(readSeq, basePath[currNodeId].Seq...)
			}
		}

		var writePath = make([]uint32, 0)
		for j = 0; j < len(readPath); j++ {
			writePath = append(writePath, readPath[j].Id)
		}
		current = &sam.SamAln{}
		current.RName = readPath[0].Name
		current.Seq = readSeq
		current.Pos = int64(readStartPos) + 1
		current.Extra = simpleGraph.PathToString(writePath, graph)
		current.Cigar = []*cigar.Cigar{&cigar.Cigar{int64(readLen), 'M'}} // TODO: trouble on this line
		answer = append(answer, current)
	}
	return answer
}

func TestGraphVariants(t *testing.T) {
	answer := GraphVariants(MakeTestGraph(), "testdata2", 0, 1, 0, 1, 0, 10, false)

	if answer == nil || len(answer) == 0 {
		t.Errorf("Problem with GraphVariants")
	}

	vcf.PrintVcf(answer)
}
