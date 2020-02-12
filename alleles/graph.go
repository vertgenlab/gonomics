package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"io/ioutil"
	"log"
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
		current.Pos = int64(readStartPos) + 1
		current.Extra = simpleGraph.PathToString(writePath, graph)
		current.Cigar = []*cigar.Cigar{&cigar.Cigar{int64(readLen), 'M'}} // TODO: trouble on this line
		answer = append(answer, current)
	}
	return answer
}



// Functions to find vars that exist in graph
type SampleData struct {
	Sample	string
	Counts 	[]uint32
}

type BatchData struct {
	Samples []string
	Counts 	[][]uint32
}

type NodeP struct {
	Sample 	string
	Node 	*simpleGraph.Node
	pVal	float64
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

func addPathToCount(data []uint32, path []uint32) {
	for i := 0; i < len(path); i++ {
		data[path[i]]++
	}
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

func CountPathsInDir(graph *simpleGraph.SimpleGraph, inDirectory string) *BatchData {
	var wg sync.WaitGroup
	receivePathCount := make(chan *SampleData)

	if !strings.HasSuffix(inDirectory, "/") {
		inDirectory = inDirectory + "/"
	}

	files, _ := ioutil.ReadDir(inDirectory)

	for _, file := range files {
		filePath := fmt.Sprintf("%s%s", inDirectory, file.Name())
		// Count Paths for Each Sample in Directory
		wg.Add(1)
		go func(filePath string, fileName string) {
			pathCount := CountPaths(graph, filePath)
			tmp := &SampleData{fileName, pathCount}
			receivePathCount <- tmp
			wg.Done()
		}(filePath, file.Name())
	}

	go func() {
		wg.Wait()
		close(receivePathCount)
	}()

	// Get path counts for each sample
	pathCounts := make([][]uint32, 0)
	samples := make([]string, 0)
	for count := range receivePathCount {
		samples = append(samples, count.Sample)
		pathCounts = append(pathCounts, count.Counts)
	}

	answer := &BatchData{samples, pathCounts}
	return answer
}

func findMajorLocalNode(nodesToCheck []*simpleGraph.Node, pathCounts [][]uint32) *simpleGraph.Node {
	var answer *simpleGraph.Node

	var i, j, sumCount, maxCount int

	maxCount = 0
	for i = 0; i < len(nodesToCheck); i++ {
		sumCount = 0
		for j = 0; j < len(pathCounts); j++ {
			sumCount += int(pathCounts[j][nodesToCheck[i].Id])
		}
		if sumCount > maxCount {
			maxCount = sumCount
			answer = nodesToCheck[i]
		}
	}

	if maxCount == 0 {
		return nil
	} else {
		return answer
	}
}

func FindMajorLocalNode(start *simpleGraph.Node, pathCounts [][]uint32) (*simpleGraph.Node, *simpleGraph.Edge) {
	var i int
	var maxNode, maxPrevNode *simpleGraph.Node
	var maxEdge *simpleGraph.Edge

	var prevNodes, nextNodes []*simpleGraph.Node

	for i = 0; i < len(start.Prev); i++ {
		prevNodes = append(prevNodes, start.Prev[i].Dest)
	}

	for i = 0; i < len(start.Next); i++ {
		nextNodes = append(nextNodes, start.Next[i].Dest)
	}

	maxNode = findMajorLocalNode(append(prevNodes, nextNodes...), pathCounts)
	maxPrevNode = findMajorLocalNode(prevNodes, pathCounts)


	for i = 0; i < len(start.Prev); i++ {
		if start.Prev[i].Dest == maxPrevNode {
			maxEdge = start.Prev[i]
			break
		}
	}

	if maxNode == nil {
		return nil, nil
	} else {
		return maxNode, maxEdge
	}
}

func calcRareNode(wg *sync.WaitGroup, sendResult chan *NodeP, start *simpleGraph.Node, sampleData *BatchData, maxPopFreq float64, minReadFreq float64, minPval float64) {
	defer wg.Done()
	pathCounts := sampleData.Counts

	majorLocalNode, majorPrevEdge := FindMajorLocalNode(start, pathCounts)

	// If there were no reads in any adjacent nodes then exit
	if majorLocalNode == nil {
		return
	}

	// If the prob of entry from the majorLocalNode is > maxPopFreq then exit.
	// No probabilty for first node in graph


	var prob float64
	if len(start.Prev) > 0 {
		prob = majorPrevEdge.Prob
	} else {
		prob = 1
	}

	if prob > maxPopFreq {
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
			b += int(pathCounts[j][majorLocalNode.Id])
			d += int(pathCounts[j][start.Id])
		}
		b = b - a
		d = d - c

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
			answer := &NodeP{sampleData.Samples[i], start, p}
			sendResult <- answer
		}
	}
}

func FindRareNodes(graph *simpleGraph.SimpleGraph, sampleData *BatchData, maxPopFreq float64, minReadFreq float64, minPval float64) []*NodeP {
	var answer = make([]*NodeP, 0)
	result := make(chan *NodeP)
	var wg sync.WaitGroup
	var i, j int

	pathCounts := sampleData.Counts

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
		go calcRareNode(&wg, result, graph.Nodes[i], sampleData, maxPopFreq, minReadFreq, minPval)
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
		Qual: 	0,
		Filter: "",
		Info: 	fmt.Sprintf("Sample:%s, p:%e", node.Sample, node.pVal),
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
		currNode := nodes[i]
		go func() {
			receiveVcf <- nodeToVcf(currNode)
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

func FindVarsInsideGraph(graph *simpleGraph.SimpleGraph, sampleData *BatchData, maxPopFreq float64, minReadFreq float64, minPval float64) []*vcf.Vcf {
	var answer []*vcf.Vcf

	// Find rare nodes in graph
	rareNodes := FindRareNodes(graph, sampleData, maxPopFreq, minReadFreq, minPval)

	answer = NodesToVcf(rareNodes)

	return answer
}


// Functions to find vars that DO NOT exist in graph
type GraphSampleMap map[GraphLocation]*AlleleCount

type BatchGraphSampleMap map[GraphLocation][]*BatchAlleleCount

type GraphLocation struct {
	Node 	*simpleGraph.Node
	Pos 	int64
}

func GraphCountAlleles(graph *simpleGraph.SimpleGraph, samFilename string, minMapQ int64) GraphSampleMap {
	AlleleMap := make(GraphSampleMap)

	var i, k int32
	var j int
	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	var done = false
	var RefIndex, SeqIndex int64
	var currentSeq []dna.Base
	var aln *sam.SamAln
	var currentIndel Indel
	var indelSeq []dna.Base
	var OrigRefIndex int64
	var Match bool

	sam.ReadHeader(samFile)

	var progressMeter int32
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {

		readPath := StringToPath(aln.Extra)

		if progressMeter%10000 == 0 {
			log.Printf("#Read %d Alignments\n", progressMeter)
		}
		progressMeter++

		// If read is unmapped then go to the next alignment
		if aln.Cigar[0].Op == '*' {
			continue
		}

		// If mapping quality is less than the threshold then go to next alignment
		if aln.MapQ < minMapQ {
			continue
		}

		SeqIndex = 0
		RefIndex = aln.Pos - 1

		currNode := 0
		ref := graph.Nodes[readPath[currNode]]

		for i = 0; i < int32(len(aln.Cigar)); i++ {
			currentSeq = aln.Seq

			//Handle deletion relative to ref
			//Each position deleted is annotated with counts + 1
			if aln.Cigar[i].Op == 'D' {
				OrigRefIndex = RefIndex
				OrigNode := ref
				indelSeq = make([]dna.Base, 1)

				// First base in indel is the base prior to the indel sequence per VCF standard format
				indelSeq[0] = OrigNode.Seq[OrigRefIndex-1]

				for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

					// If the position has already been added to the map, move along
					_, ok := AlleleMap[GraphLocation{ref, RefIndex}]

					// If the position is NOT in the map, initialize
					if !ok {
						AlleleMap[GraphLocation{ref, RefIndex}] = &AlleleCount{
							Ref: ref.Seq[RefIndex], Counts: 0, BaseA: make([]int32, 3), BaseC: make([]int32, 3), BaseG: make([]int32, 3), BaseT: make([]int32, 3), Indel: make([]Indel, 0)}
					}

					// Keep track of deleted sequence
					indelSeq = append(indelSeq, ref.Seq[RefIndex])

					AlleleMap[GraphLocation{ref, RefIndex}].Counts++

					if RefIndex + 1 == int64(len(ref.Seq)) {
						RefIndex = 0
						currNode++
						if currNode + 1 <= len(readPath) {
							ref = graph.Nodes[readPath[currNode]]
						} else {break}
					} else {
						RefIndex++
					}

				}

				Match = false
				for j = 0; j < len(AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel); j++ {
					// If the deletion has already been seen before, increment the existing entry
					// For a deletion the indelSeq should match the Ref
					if dna.CompareSeqsIgnoreCase(indelSeq, AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel[j].Ref) == 0 &&
						dna.CompareSeqsIgnoreCase(indelSeq[:1], AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel[j].Alt) == 0 {
						AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel[j].Count[0]++
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel[j].Count[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel[j].Count[2]++
						}

						Match = true
						break
					}
				}

				// If the deletion has not been seen before, then append it to the Del slice
				// For Alt indelSeq[:1] is used to give me a slice of just the first base in the slice which we defined earlier
				if Match == false {

					currentIndel = Indel{indelSeq, indelSeq[:1], make([]int32, 3)}
					currentIndel.Count[0]++
					if sam.IsForwardRead(aln) == true {
						currentIndel.Count[1]++
					} else if sam.IsReverseRead(aln) == false {
						currentIndel.Count[2]++
					}
					AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel = append(AlleleMap[GraphLocation{OrigNode, OrigRefIndex}].Indel, currentIndel)
				}

				//Handle insertion relative to ref
				//The base after the inserted sequence is annotated with an Ins read
			} else if aln.Cigar[i].Op == 'I' {

				// If the position has already been added to the map, move along
				_, ok := AlleleMap[GraphLocation{ref, RefIndex}]

				// If the position is NOT in the map, initialize
				if !ok {
					AlleleMap[GraphLocation{ref, RefIndex}] = &AlleleCount{
						Ref: ref.Seq[RefIndex], Counts: 0, BaseA: make([]int32, 3), BaseC: make([]int32, 3), BaseG: make([]int32, 3), BaseT: make([]int32, 3), Indel: make([]Indel, 0)}
				}

				// Loop through read sequence and keep track of the inserted bases
				indelSeq = make([]dna.Base, 1)

				// First base in indel is the base prior to the indel sequence per VCF standard format
				indelSeq[0] = ref.Seq[RefIndex-1]

				// Keep track of inserted sequence by moving along the read
				for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {
					indelSeq = append(indelSeq, currentSeq[SeqIndex])
					SeqIndex++
				}

				Match = false
				for j = 0; j < len(AlleleMap[GraphLocation{ref, RefIndex}].Indel); j++ {
					// If the inserted sequence matches a previously inserted sequence, then increment the count
					// For an insertion, the indelSeq should match the Alt
					if dna.CompareSeqsIgnoreCase(indelSeq, AlleleMap[GraphLocation{ref, RefIndex}].Indel[j].Alt) == 0 &&
						dna.CompareSeqsIgnoreCase(indelSeq[:1], AlleleMap[GraphLocation{ref, RefIndex}].Indel[j].Ref) == 0 {
						AlleleMap[GraphLocation{ref, RefIndex}].Indel[j].Count[0]++
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].Indel[j].Count[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].Indel[j].Count[2]++
						}
						Match = true
						break
					}
				}

				if Match == false {
					currentIndel = Indel{indelSeq[:1], indelSeq, make([]int32, 3)}
					currentIndel.Count[0]++
					if sam.IsForwardRead(aln) == true {
						currentIndel.Count[1]++
					} else if sam.IsReverseRead(aln) == true {
						currentIndel.Count[2]++
					}
					AlleleMap[GraphLocation{ref, RefIndex}].Indel = append(AlleleMap[GraphLocation{ref, RefIndex}].Indel, currentIndel)
				}

				// Note: Insertions do not contribute to the total counts as the insertion is associated with the previous reference base

				//Handle matching pos relative to ref
			} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {

				for k = 0; k < int32(aln.Cigar[i].RunLength); k++ {

					//if the position has already been added to the matrix, move along
					_, ok := AlleleMap[GraphLocation{ref, RefIndex}]

					//if the position is NOT in the matrix, add it
					if !ok {
						AlleleMap[GraphLocation{ref, RefIndex}] = &AlleleCount{
							Ref: ref.Seq[RefIndex], Counts: 0, BaseA: make([]int32, 3), BaseC: make([]int32, 3), BaseG: make([]int32, 3), BaseT: make([]int32, 3), Indel: make([]Indel, 0)}
					}

					switch currentSeq[SeqIndex] {
					case dna.A:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseA[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseA[2]++
						}
						AlleleMap[GraphLocation{ref, RefIndex}].BaseA[0]++
						AlleleMap[GraphLocation{ref, RefIndex}].Counts++
					case dna.T:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseT[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseT[2]++
						}
						AlleleMap[GraphLocation{ref, RefIndex}].BaseT[0]++
						AlleleMap[GraphLocation{ref, RefIndex}].Counts++
					case dna.G:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseG[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseG[2]++
						}
						AlleleMap[GraphLocation{ref, RefIndex}].BaseG[0]++
						AlleleMap[GraphLocation{ref, RefIndex}].Counts++
					case dna.C:
						if sam.IsForwardRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseC[1]++
						} else if sam.IsReverseRead(aln) == true {
							AlleleMap[GraphLocation{ref, RefIndex}].BaseC[2]++
						}
						AlleleMap[GraphLocation{ref, RefIndex}].BaseC[0]++
						AlleleMap[GraphLocation{ref, RefIndex}].Counts++
					}
					SeqIndex++
					if RefIndex + 1 == int64(len(ref.Seq)) {
						RefIndex = 0
						currNode++
						if currNode + 1 <= len(readPath) {
							ref = graph.Nodes[readPath[currNode]]
						} else {break}
					} else {
						RefIndex++
					}
				}
			} else if aln.Cigar[i].Op != 'H' {
				SeqIndex = SeqIndex + aln.Cigar[i].RunLength
			}
		}
	}


	return AlleleMap
}

func GraphCountAllelesInDir() []GraphSampleMap {
	answer := make([]GraphSampleMap, 0)

	return answer
}

func CreateBatchGraphSampleMap() BatchGraphSampleMap {
	answer := make(BatchGraphSampleMap)

	return answer
}

func GraphScoreVariants() {

}

func GraphScoresToVcf() []*vcf.Vcf {
	answer := make([]*vcf.Vcf, 0)

	return answer
}

func FindVarsOutsideGraph(graph *simpleGraph.SimpleGraph, simpleData *BatchData) []*vcf.Vcf {
	var answer []*vcf.Vcf

	return answer
}