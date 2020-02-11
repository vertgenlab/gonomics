package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"io/ioutil"
	"math/rand"
	"strconv"
	"strings"
	"sync"
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
func MakeTestReads(graph *simpleGraph.SimpleGraph, numReads int, readLen int, numSNP int, numSomaticSNV int) []*sam.SamAln {
	answer := make([]*sam.SamAln, 0)
	rand.Seed(time.Now().UnixNano())

	var i, j int
	var possibleSNPids = make([]int, 0)
	var SNPids = make([]int, 0)
	fmt.Println(possibleSNPids)
	// Choose random path through graph
	basePath := RandPath(graph)

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
		current.Pos = int64(readStartPos)
		current.Extra = simpleGraph.PathToString(writePath, graph)

		answer = append(answer, current)

	}
	return answer
}


// Practical Functions
type NodeP struct {
	Node 	*simpleGraph.Node
	pVal	float64
}

func addPathToCount(data []uint32, path []uint32) {
	for i := 0; i < len(path); i++ {
		data[path[i]]++
	}
}

func StringToPath(input string) []uint32 {
	answer := make([]uint32, 0)
	words := strings.Split(input, ":")
	for i := 0; i < len(words); i++ {
		node, _ := strconv.ParseUint(words[i], 10, 32)
		answer = append(answer, uint32(node))
	}
	return answer
}

func CountPaths(graph *simpleGraph.SimpleGraph, samFilename string) []uint32 {
	answer := make([]uint32, len(graph.Nodes))
	var done = false
	var aln *sam.SamAln

	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	sam.ReadHeader(samFile)

	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {
		addPathToCount(answer, StringToPath(aln.Extra))
	}
	return answer
}

func FindMajorLocalNode(start *simpleGraph.Node, pathCounts [][]uint32) (*simpleGraph.Node, *simpleGraph.Edge) {
	var i, j, sumCount, maxCount int
	var maxNode *simpleGraph.Node
	var maxEdge *simpleGraph.Edge

	var nodesToCheck = make([]*simpleGraph.Node, 0)

	for i = 0; i < len(start.Prev); i++ {
		nodesToCheck = append(nodesToCheck, start.Prev[i].Dest)
	}
	for i = 0; i < len(start.Next); i++ {
		nodesToCheck = append(nodesToCheck, start.Next[i].Dest)
	}

	maxCount = 0
	for i = 0; i < len(nodesToCheck); i++ {
		sumCount = 0
		for j = 0; j < len(pathCounts); j++ {
			sumCount += int(pathCounts[j][nodesToCheck[i].Id])
		}
		if sumCount > maxCount {
			maxNode = nodesToCheck[i]
		}
	}

	edges := append(start.Prev, start.Next...)
	for i = 0; i < len(edges); i++ {
		if edges[i].Dest.Id == start.Id {
			maxEdge = edges[i]
			break
		}
	}

	if maxCount == 0 {
		return nil, nil
	} else {
		return maxNode, maxEdge
	}
}

func calcRareNode(wg *sync.WaitGroup, sendResult chan *NodeP, start *simpleGraph.Node, pathCounts [][]uint32, maxPopFreq float64, minReadFreq float64, minPval float64) {
	defer wg.Done()

	majorLocalNode, majorLocalEdge := FindMajorLocalNode(start, pathCounts)

	// If there were no reads in any adjacent nodes then exit
	if majorLocalNode == nil {
		return
	}

	// If the prob of entry from the majorLocalNode is > maxPopFreq then exit
	if majorLocalEdge.Prob > maxPopFreq {
		return
	}

	// Begin testing on each sample
	// test is for the matrix:
	// [a b]
	// [c d]
	// a = Samples Ref Path Count
	// b = Background Ref Path Count - Samples Ref Path Count
	// c = Samples Alt Path Count
	// d = Background Alt Path Count - Samples Alt Path Count
	for i := 0; i < len(pathCounts); i++ {
		a := int(pathCounts[i][majorLocalNode.Id])
		c := int(pathCounts[i][start.Id])

		var b, d int
		for j := 0; j < len(pathCounts); j++ {
			b += int(pathCounts[j][majorLocalNode.Id]) - a
			d += int(pathCounts[j][start.Id]) - c
		}

		p := numbers.FisherExact(a, b, c, d, true)

		switch {
			// If alternate allele is zero then there is no variant and should exit
			case c == 0:
				return

			// If a = b and c = d then it is testing itself and should exit
			case a == b && c == d:
				return
		}

		if p < minPval {
			answer := &NodeP{start, p}
			sendResult <- answer
		}
	}
}

func FindRareNodes(graph *simpleGraph.SimpleGraph, pathCounts [][]uint32, maxPopFreq float64, minReadFreq float64, minPval float64) []*NodeP {
	var answer = make([]*NodeP, 0)
	result := make(chan *NodeP)
	var wg sync.WaitGroup
	var i, j int

	// Note: pathCounts is a [sample][nodes]
	// For each node ---->
	for i = 0; i < len(pathCounts[0]); i++ {

		var areThereReads = false
		// For each sample ---->
		for j = 0; j < len(pathCounts); j++ {
			if pathCounts[i][j] != 0 {
				areThereReads = true
				break
			}
		}

		// Ignore node if no reads are mapped to it
		if areThereReads == false {
			continue
		}


		wg.Add(1)
		go calcRareNode(&wg, result, graph.Nodes[i], pathCounts, maxPopFreq, minReadFreq, minPval)
	}

	go func() {
		wg.Wait()
		close(result)
	}()

	for node := range result {
		answer = append(answer, node)
	}
	return answer
}

func nodeToVcf(node *NodeP) *vcf.Vcf {
	answer := &vcf.Vcf{
		Chr: 	node.Node.Name,
		Pos: 	0,
		Id: 	"",
		Ref: 	"",
		Alt: 	dna.BasesToString(node.Node.Seq),
		Qual: 	node.pVal,
		Filter: "",
		Info: 	"",
		Format: "",
		Sample: []string{}}

	return answer
}

func NodesToVcf(nodes []*NodeP) []*vcf.Vcf {
	answer := make([]*vcf.Vcf, 0)
	var wg sync.WaitGroup
	var receiveVcf = make(chan *vcf.Vcf)

	for i := 0; i < len(nodes); i++ {
		wg.Add(1)
		go func() {
			receiveVcf <- nodeToVcf(nodes[i])
			wg.Done()
		}()
	}

	go func() {
		wg.Wait()
		close(receiveVcf)
	}()

	for record := range receiveVcf {
		answer = append(answer, record)
	}

	return answer
}

func FindVarsInGraph(graph *simpleGraph.SimpleGraph, inDirectory string, maxPopFreq float64, minReadFreq float64, minPval float64) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var wg sync.WaitGroup
	receivePathCount := make(chan []uint32)

	if !strings.HasSuffix(inDirectory, "/") {
		inDirectory = inDirectory + "/"
	}

	files, _ := ioutil.ReadDir(inDirectory)

	for _, file := range files {
		// Count Paths for Each Sample in Directory
		wg.Add(1)
		go func() {
			pathCount := CountPaths(graph, file.Name())
			receivePathCount <- pathCount
			wg.Done()
		}()
	}

	go func() {
		wg.Wait()
		close(receivePathCount)
	}()

	// Get path counts for each sample
	pathCounts := make([][]uint32, 0)
	for count := range receivePathCount {
		pathCounts = append(pathCounts, count)
	}

	// Find rare nodes in graph
	rareNodes := FindRareNodes(graph, pathCounts, maxPopFreq, minReadFreq, minPval)

	answer = NodesToVcf(rareNodes)

	return answer
}