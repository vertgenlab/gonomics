package genomeGraph

import (
	"bytes"
	"io"
	"log"
	"math"
	"sort"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers"
)

const (
	defaultMatrixSize int  = 2480
	LeftDirection     byte = 0
	RightDirection    byte = 1
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
}

type dynamicScoreKeeper struct {
	i       int
	j       int
	currMax int64
	route   []cigar.ByteCigar
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

func resetDynamicScore(sk dynamicScoreKeeper) {
	sk.route = sk.route[:0]
	sk.currMax = 0
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

func resetScoreKeeper(sk scoreKeeper) {
	sk.targetStart, sk.targetEnd = 0, 0
	sk.queryStart, sk.queryEnd = 0, 0
	sk.currScore, sk.seedScore = 0, 0
	sk.perfectScore = 0
	sk.leftAlignment, sk.rightAlignment = sk.leftAlignment[:0], sk.rightAlignment[:0]
	sk.leftPath, sk.rightPath = sk.leftPath[:0], sk.rightPath[:0]
	sk.leftSeq, sk.rightSeq, sk.currSeq = sk.leftSeq[:0], sk.rightSeq[:0], sk.currSeq[:0]
	sk.leftScore, sk.rightScore = 0, 0
}

// AlignGraphTraversal performs graph traversal for alignment either to the left or right, calculating optimal alignment using dynamic programming.
func AlignGraphTraversal(n *Node, seq []dna.Base, position int, currentPath []uint32, extension int, read []dna.Base, scores [][]int64, matrix *MatrixAln, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool, direction byte) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	s := pool.Get().(*dnaPool)
	defer pool.Put(s)

	s.Seq, s.Path = s.Seq[:0], s.Path[:0]
	s.Path = append([]uint32(nil), currentPath...)
	if direction == LeftDirection {
		s.Seq = nodeSeqTraversal(n, extension, position, seq, s.Seq, LeftDirection)
		AddPath(s.Path, n.Id)
	} else {
		s.Seq = nodeSeqTraversal(n, extension, position, seq, s.Seq, RightDirection)
	}

	var alignment []cigar.ByteCigar
	var score int64
	var targetStart, queryStart int
	if direction == LeftDirection {
		if len(s.Seq) >= extension || len(n.Prev) == 0 {
			score, alignment, targetStart, queryStart = DynamicAln(s.Seq, read, scores, matrix, -600, dynamicScore, LeftDirection)
			targetStart = position - len(s.Seq) + targetStart
		} else {
			score = math.MinInt64
			for _, prevNode := range n.Prev {
				route, currScore, targetStart, queryStart, path := AlignGraphTraversal(prevNode.Dest, s.Seq, len(prevNode.Dest.Seq), s.Path, extension, read, scores, matrix, sk, dynamicScore, pool, LeftDirection)
				if currScore > score {
					score, alignment, sk.targetStart, sk.queryStart, sk.leftPath = currScore, route, targetStart, queryStart, path
				}
			}
			cigar.ReverseBytesCigar(alignment)
			ReversePath(sk.leftPath)
		}
	} else {
		if len(s.Seq) >= extension || len(n.Next) == 0 {
			score, alignment, targetStart, queryStart = DynamicAln(s.Seq, read, scores, matrix, -600, dynamicScore, RightDirection)
			targetStart += position
		} else {
			score = math.MinInt64
			for _, nextNode := range n.Next {
				route, currScore, targetEnd, queryEnd, path := AlignGraphTraversal(nextNode.Dest, s.Seq, 0, s.Path, extension, read, scores, matrix, sk, dynamicScore, pool, RightDirection)
				if currScore > score {
					score, alignment, sk.targetEnd, sk.queryEnd, sk.rightPath = currScore, route, targetEnd, queryEnd, path
				}
			}
		}
	}

	return alignment, score, targetStart, queryStart, s.Path
}

// DynamicAln performs dynamic programming alignment with traceback, adaptable for both leftward and rightward directions.
func DynamicAln(alpha, beta []dna.Base, scores [][]int64, matrix *MatrixAln, gapPen int64, dynamicScore dynamicScoreKeeper, direction byte) (int64, []cigar.ByteCigar, int, int) {
	resetDynamicScore(dynamicScore)
	var i, j int
	var rows, columns = len(alpha), len(beta)

	// Initialize the matrix's first row and column based on gap penalties.
	for i = 0; i <= rows; i++ {
		matrix.m[i][0] = int64(i) * gapPen
		matrix.trace[i][0] = 'D' // Deletion
	}
	for j = 1; j <= columns; j++ {
		matrix.m[0][j] = int64(j) * gapPen
		matrix.trace[0][j] = 'I' // Insertion
	}

	// Fill the matrix using dynamic programming.
	for i = 1; i <= rows; i++ {
		for j = 1; j <= columns; j++ {
			matchScore := matrix.m[i-1][j-1] + scores[alpha[i-1]][beta[j-1]]
			insertScore := matrix.m[i][j-1] + gapPen
			deleteScore := matrix.m[i-1][j] + gapPen

			matrix.m[i][j], matrix.trace[i][j] = maxScore(matchScore, insertScore, deleteScore)

			// Update the maximum score and its position.
			if matrix.m[i][j] > dynamicScore.currMax {
				dynamicScore.currMax = matrix.m[i][j]
				dynamicScore.i, dynamicScore.j = i, j
			}
		}
	}

	// Traceback to construct the alignment.
	alignment, alignStartI, alignStartJ := traceback(matrix, dynamicScore, direction)
	return matrix.m[dynamicScore.i][dynamicScore.j], alignment, alignStartI, alignStartJ
}

// maxScore selects the maximum score and corresponding trace for dynamic programming.
func maxScore(matchScore, insertScore, deleteScore int64) (int64, byte) {
	maxScore := matchScore
	trace := cigar.Match
	if insertScore > maxScore {
		maxScore = insertScore
		trace = cigar.Insertion // Insertion
	}
	if deleteScore > maxScore {
		maxScore = deleteScore
		trace = cigar.Deletion // Deletion
	}
	return maxScore, trace
}

// traceback constructs the alignment from the dynamic programming matrix, starting from the maximum score position.
func traceback(matrix *MatrixAln, dynamicScore dynamicScoreKeeper, direction byte) ([]cigar.ByteCigar, int, int) {
	var alignment []cigar.ByteCigar
	i, j := dynamicScore.i, dynamicScore.j

	// The traceback process might be adjusted based on the alignment direction.
	// For the sake of this example, assume it's the same for both directions.
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

	// For leftward alignment, additional adjustments might be needed.
	if direction == LeftDirection {
		// Adjust alignment start indices if needed.
		i, j = adjustAlignmentStart(matrix, i, j)

	}

	return alignment, i, j // Return the starting indices of the alignment for further processing.
}

func nodeSeqTraversal(n *Node, extension int, position int, seq, ans []dna.Base, direction byte) []dna.Base {
	switch direction {
	case LeftDirection:
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
	case RightDirection:
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

// adjustAlignmentStart finds the starting indices of the alignment by tracing back from the given indices until a non-positive score is encountered
func adjustAlignmentStart(matrix *MatrixAln, i, j int) (int, int) {
	for i > 0 && j > 0 && matrix.m[i][j] > 0 {
		switch matrix.trace[i][j] {
		case cigar.Match, cigar.Insertion, cigar.Deletion:
			i--
			j--
		}
	}
	return i, j
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

func extendToTheLeftDev(node *Node, read *fastq.FastqBig, currPart SeedDev) []SeedDev {
	var answer, prevParts []SeedDev
	var i int
	var readBase dna.Base

	if currPart.QueryStart > 0 && currPart.TargetStart == 0 {
		for i = 0; i < len(node.Prev); i++ {
			if currPart.PosStrand {
				readBase = dnaTwoBit.GetBase(&read.Rainbow[0], uint(currPart.QueryStart)-1)
			} else {
				readBase = dnaTwoBit.GetBase(&read.RainbowRc[0], uint(currPart.QueryStart)-1)
			}
			if readBase == dnaTwoBit.GetBase(node.Prev[i].Dest.SeqTwoBit, uint(node.Prev[i].Dest.SeqTwoBit.Len)-1) {
				prevParts = extendToTheLeftHelperDev(node.Prev[i].Dest, read, currPart)
				answer = append(answer, prevParts...)
			}
		}
	}

	if len(answer) == 0 {
		return []SeedDev{currPart}
	} else {
		return answer
	}
}

func extendToTheLeftHelperDev(node *Node, read *fastq.FastqBig, nextPart SeedDev) []SeedDev {
	const basesPerInt int = 32
	var nodePos int = node.SeqTwoBit.Len - 1
	var readPos int = int(nextPart.QueryStart) - 1
	var nodeOffset int = nodePos % basesPerInt
	var readOffset int = 31 - ((readPos - nodeOffset + 31) % 32)
	var leftMatches int = 0
	var currPart SeedDev
	var prevParts, answer []SeedDev
	var i int
	var readBase dna.Base

	if nextPart.PosStrand {
		leftMatches = numbers.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.Rainbow[readOffset], readPos+readOffset))
	} else {
		leftMatches = numbers.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.RainbowRc[readOffset], readPos+readOffset))
	}

	if leftMatches == 0 {
		log.Fatal("Error: should not have zero matches to the left\n")
	}
	currPart = SeedDev{TargetId: node.Id, TargetStart: uint32(nodePos - (leftMatches - 1)), QueryStart: uint32(readPos - (leftMatches - 1)), Length: uint32(leftMatches), PosStrand: nextPart.PosStrand, TotalLength: uint32(leftMatches) + nextPart.TotalLength, NextPart: &nextPart}
	// we went all the way to end and there might be more
	if currPart.QueryStart > 0 && currPart.TargetStart == 0 {
		for i = 0; i < len(node.Prev); i++ {
			if nextPart.PosStrand {
				readBase = dnaTwoBit.GetBase(&read.Rainbow[0], uint(currPart.QueryStart)-1)
			} else {
				readBase = dnaTwoBit.GetBase(&read.RainbowRc[0], uint(currPart.QueryStart)-1)
			}
			if readBase == dnaTwoBit.GetBase(node.Prev[i].Dest.SeqTwoBit, uint(node.Prev[i].Dest.SeqTwoBit.Len)-1) {
				prevParts = extendToTheLeftHelperDev(node.Prev[i].Dest, read, currPart)
				answer = append(answer, prevParts...)
			}
		}
	}
	// if the alignment did not go to another node, return the match for this node
	if len(answer) == 0 {
		answer = []SeedDev{currPart}
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

func restartSeedHelper(helper *seedHelper) {
	helper.currHits = helper.currHits[:0]
	helper.keyIdx, helper.keyOffset, helper.readOffset, helper.nodeOffset = 0, 0, 0, 0
	helper.nodeIdx, helper.nodePos = 0, 0
	helper.seqKey, helper.codedNodeCoord = 0, 0
	helper.leftMatches = 0
}

// seedBuildHelper.nodeIdx, seedBuildHelper.nodePos int64 = 0, 0.
func seedMapMemPool(seedHash map[uint64][]uint64, nodes []Node, read *fastq.FastqBig, seedLen int, perfectScore int64, scoreMatrix [][]int64, finalSeeds []SeedDev, tempSeeds []SeedDev, seedBuildHelper *seedHelper) []SeedDev {
	const basesPerInt int64 = 32
	restartSeedHelper(seedBuildHelper)
	seedBuildHelper.keyShift = 64 - (uint(seedLen) * 2)

	for readStart := 0; readStart < len(read.Seq)-seedLen+1; readStart++ {
		seedBuildHelper.keyIdx = (readStart + 31) / 32
		seedBuildHelper.keyOffset = 31 - ((readStart + 31) % 32)
		// do fwd strand
		seedBuildHelper.seqKey = read.Rainbow[seedBuildHelper.keyOffset].Seq[seedBuildHelper.keyIdx] >> seedBuildHelper.keyShift
		seedBuildHelper.currHits = seedHash[seedBuildHelper.seqKey]

		for _, seedBuildHelper.codedNodeCoord = range seedBuildHelper.currHits {
			seedBuildHelper.nodeIdx, seedBuildHelper.nodePos = numberToChromAndPos(seedBuildHelper.codedNodeCoord)
			seedBuildHelper.nodeOffset = int(seedBuildHelper.nodePos % basesPerInt)
			seedBuildHelper.readOffset = 31 - ((readStart - seedBuildHelper.nodeOffset + 31) % 32)
			seedBuildHelper.leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[seedBuildHelper.nodeIdx].SeqTwoBit, int(seedBuildHelper.nodePos), &read.Rainbow[seedBuildHelper.readOffset], readStart+seedBuildHelper.readOffset))
			tempSeeds = extendToTheRightDev(&nodes[seedBuildHelper.nodeIdx], read, readStart-(seedBuildHelper.leftMatches-1), int(seedBuildHelper.nodePos)-(seedBuildHelper.leftMatches-1), true, tempSeeds)
			for _, seedBuildHelper.tempSeed = range tempSeeds {
				finalSeeds = append(finalSeeds, extendToTheLeftDev(&nodes[seedBuildHelper.nodeIdx], read, seedBuildHelper.tempSeed)...)
			}
		}
		// do rev strand
		seedBuildHelper.seqKey = read.RainbowRc[seedBuildHelper.keyOffset].Seq[seedBuildHelper.keyIdx] >> seedBuildHelper.keyShift
		seedBuildHelper.currHits = seedHash[seedBuildHelper.seqKey]
		for _, seedBuildHelper.codedNodeCoord = range seedBuildHelper.currHits {
			seedBuildHelper.nodeIdx, seedBuildHelper.nodePos = numberToChromAndPos(seedBuildHelper.codedNodeCoord)
			seedBuildHelper.nodeOffset = int(seedBuildHelper.nodePos % basesPerInt)
			seedBuildHelper.readOffset = 31 - ((readStart - seedBuildHelper.nodeOffset + 31) % 32)

			seedBuildHelper.leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[seedBuildHelper.nodeIdx].SeqTwoBit, int(seedBuildHelper.nodePos), &read.RainbowRc[seedBuildHelper.readOffset], readStart+seedBuildHelper.readOffset))
			tempSeeds = extendToTheRightDev(&nodes[seedBuildHelper.nodeIdx], read, readStart-(seedBuildHelper.leftMatches-1), int(seedBuildHelper.nodePos)-(seedBuildHelper.leftMatches-1), false, tempSeeds)
			finalSeeds = append(finalSeeds, tempSeeds...)
		}
	}
	if len(finalSeeds) > 100 {
		SortSeedLen(finalSeeds)
	} else {
		heapSortSeeds(finalSeeds)
	}
	return finalSeeds
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

func SimpleWriteGirafPair(filename string, input <-chan giraf.GirafPair, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	var buf *bytes.Buffer
	var simplePool = sync.Pool{
		New: func() interface{} {
			return &bytes.Buffer{}
		},
	}
	var err error
	for gp := range input {
		buf = simplePool.Get().(*bytes.Buffer)
		buf.Reset()
		_, err = buf.WriteString(giraf.ToString(&gp.Fwd))
		exception.PanicOnErr(err)
		err = buf.WriteByte('\n')
		exception.PanicOnErr(err)
		_, err = buf.WriteString(giraf.ToString(&gp.Rev))
		exception.PanicOnErr(err)
		err = buf.WriteByte('\n')
		exception.PanicOnErr(err)
		_, err = io.Copy(file, buf)
		exception.PanicOnErr(err)
		simplePool.Put(buf)
	}
	err = file.Close()
	exception.PanicOnErr(err)
	wg.Done()
}
