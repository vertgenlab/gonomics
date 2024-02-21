package genomeGraph

import (
	"math"
	"sort"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/numbers"
)

const (
	defaultMatrixSize int  = 2480
	leftTraversal     byte = 0
	rightTraversal    byte = 1
)

type memoryPool struct {
	Hits   []SeedDev
	Worker []SeedDev
}

type MatrixAln struct {
	m     [][]int64
	trace [][]byte
}

type dnaPool struct {
	Seq         []dna.Base
	Path        []uint32
	queryStart  int
	queryEnd    int
	targetStart int
	targetEnd   int

	currScore int64
}

type dynamicScoreKeeper struct {
	i       int
	j       int
	currMax int64
}

type scoreKeeper struct {
	targetStart  int
	targetEnd    int
	queryStart   int
	queryEnd     int
	extension    int
	currScore    int64
	seedScore    int64
	perfectScore int64
	leftScore    int64
	rightScore   int64
	leftPath     []uint32
	rightPath    []uint32
	leftSeq      []dna.Base
	rightSeq     []dna.Base
	currSeq      []dna.Base
	tailSeed     SeedDev

	currSeed       SeedDev
	leftAlignment  []cigar.ByteCigar
	rightAlignment []cigar.ByteCigar
}

type seedHelper struct {
	currHits                                  []uint64
	codedNodeCoord                            uint64
	seqKey                                    uint64
	keyShift                                  uint
	keyIdx, keyOffset, readOffset, nodeOffset int
	nodeIdx, nodePos                          int64
	leftMatches                               int
	tempSeed                                  SeedDev
}

// AlignGraphTraversal performs graph traversal for alignment either to the left or right, calculating optimal alignment using dynamic programming.
func AlignGraphTraversal(n *Node, seq []dna.Base, position int, currentPath []uint32, extension int, read []dna.Base, scores [][]int64, matrix *MatrixAln, sk *scoreKeeper, dynamicScore *dynamicScoreKeeper, pool *sync.Pool, direction byte) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	// Fetch a pooled object to reduce allocations.
	s := pool.Get().(*dnaPool)
	defer pool.Put(s)

	// Initialize sequence and path slices without reallocation.
	s.Seq = append(s.Seq[:0], seq...)
	s.Path = append(s.Path[:0], currentPath...)

	// Traverse node sequence and update the path accordingly.
	s.Seq = nodeSeqTraversal(n, extension, position, seq, s.Seq, direction)
	if direction == leftTraversal {
		s.Path = append(s.Path, n.Id) // Add node ID for left direction.
	}

	// Initialize local variables to hold alignment details.
	var alignment []cigar.ByteCigar
	var score int64
	var targetStart, queryStart int

	// Check termination condition for recursion.
	isEndOfTraversal := (direction == leftTraversal && (len(s.Seq) >= extension || len(n.Prev) == 0)) ||
		(direction == rightTraversal && (len(s.Seq) >= extension || len(n.Next) == 0))

	if isEndOfTraversal {
		// Perform dynamic alignment at the end of traversal.
		score, alignment, targetStart, queryStart = DynamicAln(s.Seq, read, scores, matrix, -600, dynamicScore, direction)
		if direction == rightTraversal {
			targetStart += position // Adjust target start for right direction.
		}
		return alignment, score, targetStart, queryStart, append([]uint32(nil), s.Path...)
	}

	// Recursive case: traverse next or previous nodes based on direction.
	score = math.MinInt64
	nodes := n.Next
	if direction == leftTraversal {
		nodes = n.Prev
	}
	for _, edge := range nodes {
		route, currScore, end, start, path := AlignGraphTraversal(edge.Dest, s.Seq, 0, s.Path, extension, read, scores, matrix, sk, dynamicScore, pool, direction)
		if currScore > score {
			// Update alignment details if current score is better.
			score, alignment, targetStart, queryStart = currScore, route, end, start
			s.Path = path
		}
	}

	// Adjust the target start for left direction after all recursive calls.
	if direction == leftTraversal {
		targetStart = position - len(s.Seq) + targetStart
		cigar.ReverseBytesCigar(alignment)
		ReversePath(s.Path)
	} else {
		targetStart += position // Adjust for right direction.
	}

	return alignment, score, targetStart, queryStart, append([]uint32(nil), s.Path...)
}

// DynamicAln performs dynamic programming-based sequence alignment by traversing a graph of nodes
func DynamicAln(alpha, beta []dna.Base, scores [][]int64, matrix *MatrixAln, gapPen int64, dynamicScore *dynamicScoreKeeper, direction byte) (int64, []cigar.ByteCigar, int, int) {
	dynamicScore.currMax = 0 // Reset the current max score.
	var i, j int
	var rows, columns int = len(alpha), len(beta)
	// Initialize the first row and column based on gap penalties.
	for i = 1; i <= rows; i++ {
		matrix.m[i][0] = int64(i) * gapPen
		matrix.trace[i][0] = 'D' // Indicate a deletion.
	}
	for j = 1; j <= columns; j++ {
		matrix.m[0][j] = int64(j) * gapPen
		matrix.trace[0][j] = 'I' // Indicate an insertion.
	}

	// Fill in the scoring and traceback matrices.
	for i = 1; i <= rows; i++ {
		for j = 1; j <= columns; j++ {
			matchOrMismatch := matrix.m[i-1][j-1] + scores[alpha[i-1]][beta[j-1]]
			deletion := matrix.m[i-1][j] + gapPen
			insertion := matrix.m[i][j-1] + gapPen

			// Choose the best score.
			if matchOrMismatch >= deletion && matchOrMismatch >= insertion {
				matrix.m[i][j] = matchOrMismatch
				matrix.trace[i][j] = cigar.Match
			} else if deletion > insertion {
				matrix.m[i][j] = deletion
				matrix.trace[i][j] = cigar.Deletion
			} else {
				matrix.m[i][j] = insertion
				matrix.trace[i][j] = cigar.Insertion
			}
			// Update the current maximum score.
			if matrix.m[i][j] > dynamicScore.currMax {
				dynamicScore.currMax = matrix.m[i][j]
				dynamicScore.i, dynamicScore.j = i, j
			}
		}
	}

	// Traceback from the maximum score's position to build the alignment.
	alignment := make([]cigar.ByteCigar, 0)
	i, j = dynamicScore.i, dynamicScore.j
	for i > 0 && j > 0 {
		switch matrix.trace[i][j] {
		case cigar.Match:
			alignment = cigar.AddCigarByte(alignment, cigar.ByteCigar{RunLen: 1, Op: cigar.Match})
			i--
			j--
		case cigar.Insertion:
			alignment = cigar.AddCigarByte(alignment, cigar.ByteCigar{RunLen: 1, Op: cigar.Insertion})
			j--
		case cigar.Deletion:
			alignment = cigar.AddCigarByte(alignment, cigar.ByteCigar{RunLen: 1, Op: cigar.Deletion})
			i--
		}
	}
	// Prepend soft clips if necessary (not shown here, depends on alignment strategy).
	return dynamicScore.currMax, alignment, dynamicScore.i, dynamicScore.j
}

func nodeSeqTraversal(n *Node, extension int, position int, seq, ans []dna.Base, direction byte) []dna.Base {
	switch direction {
	case leftTraversal:
		startPos := position - extension
		if startPos < 0 {
			startPos = 0
		}
		// Ensure the combined sequence does not exceed the original sequence length.
		endPos := position
		if endPos > len(n.Seq) {
			endPos = len(n.Seq)
		}
		// Combine sequences.
		ans = append(ans, n.Seq[startPos:endPos]...)
		ans = append(ans, seq...)
	case rightTraversal:
		endPos := position + extension
		if endPos > len(n.Seq) {
			endPos = len(n.Seq)
		}
		// Combine sequences.
		ans = append(ans, seq...)
		ans = append(ans, n.Seq[position:endPos]...)
	}
	return ans
}

func NewDnaPool() sync.Pool {
	return sync.Pool{
		New: func() interface{} {
			dnaSeq := dnaPool{
				Seq:         make([]dna.Base, 0, 150),
				Path:        make([]uint32, 0, 10),
				queryStart:  0,
				targetStart: 0,
				targetEnd:   0,
				queryEnd:    0,
			}
			return &dnaSeq
		},
	}
}

func NewMemSeedPool() sync.Pool {
	return sync.Pool{
		New: func() interface{} {
			pool := memoryPool{
				Hits:   make([]SeedDev, 0, 10000),
				Worker: make([]SeedDev, 0, 10000),
			}
			return &pool
		},
	}
}

func NewSwMatrix(size int) MatrixAln {
	sw := MatrixAln{}
	sw.m, sw.trace = MatrixSetup(size)
	return sw
}

func MatrixSetup(size int) ([][]int64, [][]byte) {
	m := make([][]int64, size)
	trace := make([][]byte, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		trace[idx] = make([]byte, size)
	}
	return m, trace
}

func ReversePath(alpha []uint32) {
	var i, off int
	for i = len(alpha)/2 - 1; i >= 0; i-- {
		off = len(alpha) - 1 - i
		alpha[i], alpha[off] = alpha[off], alpha[i]
	}
}

func seedsHeapify(a []SeedDev, i, size int) {
	largest := i
	l := 2*i + 1
	r := 2*i + 2

	if l < size && a[l].TotalLength < a[largest].TotalLength {
		largest = l
	}
	if r < size && a[r].TotalLength < a[largest].TotalLength {
		largest = r
	}

	if largest != i {
		a[i], a[largest] = a[largest], a[i]
		seedsHeapify(a, largest, size)
	}
}

func buildSeedHeap(a []SeedDev) {
	n := len(a)
	for i := n/2 - 1; i >= 0; i-- {
		seedsHeapify(a, i, n)
	}
}

func heapSortSeeds(a []SeedDev) {
	buildSeedHeap(a)
	for i := len(a) - 1; i > 0; i-- {
		a[0], a[i] = a[i], a[0]
		seedsHeapify(a, 0, i)
	}
}

func extendToTheRightDev(node *Node, read *fastq.FastqBig, readStart int, nodeStart int, posStrand bool, answer []SeedDev) []SeedDev {
	const basesPerInt int = 32
	answer = answer[:0]
	var nodeOffset int = nodeStart % basesPerInt
	var readOffset int = 31 - ((readStart - nodeOffset + 31) % 32)
	var rightMatches int = 0
	var currNode SeedDev
	var nextParts []SeedDev
	var i, j int = 0, 0
	if posStrand {
		rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, &read.Rainbow[readOffset], readStart+readOffset)
	} else {
		rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, &read.RainbowRc[readOffset], readStart+readOffset)
	}
	// nothing aligned here
	if rightMatches == 0 {
		return nil
	}
	// we went all the way to end and there might be more
	if readStart+rightMatches < len(read.Seq) && nodeStart+rightMatches == node.SeqTwoBit.Len && len(node.Next) != 0 {
		for i = 0; i < len(node.Next); i++ {
			nextParts = extendToTheRightDev(node.Next[i].Dest, read, readStart+rightMatches, 0, posStrand, nextParts)
			// if we aligned into the next node, make a seed for this node and point it to the next one
			for j = 0; j < len(nextParts); j++ {
				currNode = SeedDev{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(rightMatches), PosStrand: posStrand, TotalLength: uint32(rightMatches) + nextParts[j].TotalLength, NextPart: &nextParts[j]}
				answer = append(answer, currNode)
			}
		}
	}
	// if the alignment did not go to another node, return the match for this node
	if len(answer) == 0 {
		currNode = SeedDev{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(rightMatches), PosStrand: posStrand, TotalLength: uint32(rightMatches), NextPart: nil}
		answer = []SeedDev{currNode}
	}

	return answer
}

func newSeedBuilder() *seedHelper {
	var tmp SeedDev = SeedDev{}
	return &seedHelper{
		currHits: make([]uint64, 0, 20),
		tempSeed: tmp,
	}
}

// seedBuildHelper.nodeIdx, seedBuildHelper.nodePos int64 = 0, 0.
func seedMapMemPool(seedHash map[uint64][]uint64, nodes []Node, read *fastq.FastqBig, seedLen int,
	perfectScore int64, scoreMatrix [][]int64, finalSeeds []SeedDev, tempSeeds []SeedDev, seedBuildHelper *seedHelper) []SeedDev {

	resetSeedBuilder(seedBuildHelper, seedLen)

	for readStart := 0; readStart < len(read.Seq)-seedLen+1; readStart++ {
		getSeedKey(read, seedBuildHelper, readStart)

		processSeeds(seedHash, nodes, read, finalSeeds, tempSeeds, seedBuildHelper, readStart, true)
		processSeeds(seedHash, nodes, read, finalSeeds, tempSeeds, seedBuildHelper, readStart, false)
	}
	sortSeedsIfRequired(finalSeeds)

	return finalSeeds
}

// resetSeedBuilder resets the seedBuilder helper for use.
func resetSeedBuilder(seedBuildHelper *seedHelper, seedLen int) {
	seedBuildHelper.keyShift = 64 - (uint(seedLen) * 2)
}

// getSeedKey computes the seed key.
func getSeedKey(read *fastq.FastqBig, seedBuildHelper *seedHelper, readStart int) {
	seedBuildHelper.keyIdx = (readStart + 31) / 32
	seedBuildHelper.keyOffset = 31 - ((readStart + 31) % 32)
}

// processSeeds processes the seeds for the forward/reverse strands.
func processSeeds(seedHash map[uint64][]uint64, nodes []Node, read *fastq.FastqBig, finalSeeds []SeedDev,
	tempSeeds []SeedDev, seedBuildHelper *seedHelper, readStart int, isForward bool) {

	if isForward {
		seedBuildHelper.seqKey = read.Rainbow[seedBuildHelper.keyOffset].Seq[seedBuildHelper.keyIdx] >> seedBuildHelper.keyShift
	} else {
		seedBuildHelper.seqKey = read.RainbowRc[seedBuildHelper.keyOffset].Seq[seedBuildHelper.keyIdx] >> seedBuildHelper.keyShift
	}

	seedBuildHelper.currHits = seedHash[seedBuildHelper.seqKey]

	for _, seedBuildHelper.codedNodeCoord = range seedBuildHelper.currHits {
		processHits(nodes, read, finalSeeds, tempSeeds, seedBuildHelper, readStart, isForward)
	}
}

// processHits processes the hits for a given seed.
func processHits(nodes []Node, read *fastq.FastqBig, finalSeeds []SeedDev, tempSeeds []SeedDev, seedBuildHelper *seedHelper,
	readStart int, isForward bool) {

	seedBuildHelper.nodeIdx, seedBuildHelper.nodePos = numberToChromAndPos(seedBuildHelper.codedNodeCoord)
	seedBuildHelper.nodeOffset = int(seedBuildHelper.nodePos % basesPerInt)
	seedBuildHelper.readOffset = 31 - ((readStart - seedBuildHelper.nodeOffset + 31) % 32)

	if isForward {
		seedBuildHelper.leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[seedBuildHelper.nodeIdx].SeqTwoBit,
			int(seedBuildHelper.nodePos), &read.Rainbow[seedBuildHelper.readOffset], readStart+seedBuildHelper.readOffset))
	} else {
		seedBuildHelper.leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[seedBuildHelper.nodeIdx].SeqTwoBit,
			int(seedBuildHelper.nodePos), &read.RainbowRc[seedBuildHelper.readOffset], readStart+seedBuildHelper.readOffset))
	}

	tempSeeds = extendToTheRightDev(&nodes[seedBuildHelper.nodeIdx], read, readStart-(seedBuildHelper.leftMatches-1),
		int(seedBuildHelper.nodePos)-(seedBuildHelper.leftMatches-1), isForward, tempSeeds)

	finalSeeds = append(finalSeeds, tempSeeds...)
}

// sortSeedsIfRequired sorts the seeds if required.
func sortSeedsIfRequired(finalSeeds []SeedDev) {
	if len(finalSeeds) > 100 {
		SortSeedLen(finalSeeds)
	} else {
		heapSortSeeds(finalSeeds)
	}
}

func SortSeedLen(seeds []SeedDev) {
	sort.Slice(seeds, func(i, j int) bool { return seeds[i].TotalLength > seeds[j].TotalLength })
}

func SoftClipBases(front int, lengthOfRead int, cig []cigar.ByteCigar) []cigar.ByteCigar {
	var runLen int = cigar.QueryRunLen(cig)
	if runLen >= lengthOfRead {
		return cig
	}
	answer := make([]cigar.ByteCigar, 0, len(cig)+2)
	if front > 0 {
		answer = append(answer, cigar.ByteCigar{RunLen: uint16(front), Op: 'S'})
	}
	answer = append(answer, cig...)
	// Calculate the remaining unaligned portion at the end and add soft clipping if necessary
	remaining := lengthOfRead - (front + runLen)
	if remaining > 0 {
		answer = append(answer, cigar.ByteCigar{RunLen: uint16(remaining), Op: 'S'})
	}
	return answer
}

// func SimpleWriteGirafPair(filename string, input <-chan giraf.GirafPair, wg *sync.WaitGroup) {
// 	file := fileio.EasyCreate(filename)
// 	defer file.Close()
// 	defer wg.Done()

// 	simplePool := &sync.Pool{
// 		New: func() interface{} {
// 			return new(bytes.Buffer)
// 		},
// 	}

// 	for gp := range input {
// 		writeGirafPairToFile(gp, file, simplePool)
// 	}
// }

// func writeGirafPairToFile(gp giraf.GirafPair, file *os.File, pool *sync.Pool) {
// 	buf := pool.Get().(*bytes.Buffer)
// 	buf.Reset()

// 	writeStringAndCheckErr(buf, giraf.ToString(&gp.Fwd))
// 	writeByteAndCheckErr(buf, '\n')

// 	writeStringAndCheckErr(buf, giraf.ToString(&gp.Rev))
// 	writeByteAndCheckErr(buf, '\n')

// 	_, err := io.Copy(file, buf)
// 	exception.PanicOnErr(err)

// 	pool.Put(buf)
// }

// func WriteStringAndCheckErr(buf *bytes.Buffer, s string) {
// 	_, err := buf.WriteString(s)
// 	exception.PanicOnErr(err)
// }

// func WriteByteAndCheckErr(buf *bytes.Buffer, b byte) {
// 	err := buf.WriteByte(b)
// 	exception.PanicOnErr(err)
// }
