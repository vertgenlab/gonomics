package genomeGraph

import (
	"bytes"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/numbers"
	"io"
	"log"
	"math"
	"sort"
	"sync"
)

const (
	defaultMatrixSize int = 2480
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
	currScore   int64
}

type dynamicScoreKeeper struct {
	i        int
	j        int
	routeIdx int
	currMax  int64
	route    []cigar.ByteCigar
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
			return &dnaPool{
				Seq:  make([]dna.Base, 0, 150),
				Path: make([]uint32, 0, 10),
			}
		},
	}
}

func NewMemSeedPool() sync.Pool {
	return sync.Pool{
		New: func() interface{} {
			return &memoryPool{
				Hits:   make([]SeedDev, 0, 10000),
				Worker: make([]SeedDev, 0, 10000),
			}
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

func resetScoreKeeper(sk *scoreKeeper) {
	sk.targetStart, sk.targetEnd = 0, 0
	sk.queryStart, sk.queryEnd = 0, 0
	sk.currScore, sk.seedScore, sk.perfectScore, sk.leftScore, sk.rightScore = 0, 0, 0, 0, 0
	sk.leftAlignment = sk.leftAlignment[:0]
	sk.rightAlignment = sk.rightAlignment[:0]
	sk.leftPath = sk.leftPath[:0]
	sk.rightPath = sk.rightPath[:0]
	sk.leftSeq = sk.leftSeq[:0]
	sk.rightSeq = sk.rightSeq[:0]
	sk.currSeq = sk.currSeq[:0]
}

func getLeftTargetBases(n *Node, extension int, refEnd int, seq []dna.Base, ans []dna.Base) []dna.Base {
	toAppend := n.Seq[refEnd-numbers.Min(len(seq)+refEnd, extension) : refEnd]
	return append(ans, append(toAppend, seq...)...)
}

func getRightBases(n *Node, extension int, start int, seq []dna.Base, ans []dna.Base) []dna.Base {
	toAppend := n.Seq[start : start+numbers.Min(len(seq)+len(n.Seq)-start, extension)-len(seq)]
	return append(ans, append(seq, toAppend...)...)
}

type AlignResult struct {
	Alignment []cigar.ByteCigar
	Score     int64
	TargetPos int // It can represent either TargetStart or TargetEnd based on context.
	QueryPos  int // It can represent either QueryStart or QueryEnd based on context.
	Path      []uint32
}

func LeftAlignTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, extension int, read []dna.Base, scores [][]int64, matrix *MatrixAln, pool *sync.Pool) AlignResult {
	s := pool.Get().(*dnaPool)
	defer pool.Put(s)

	s.Seq = getLeftTargetBases(n, extension, refEnd, seq, s.Seq)
	s.Path = append(s.Path[:0], currentPath...)
	AddPath(s.Path, n.Id)

	if len(seq)+refEnd >= extension || len(n.Prev) == 0 {
		score, alignment, targetStart, queryStart := LeftDynamicAln(s.Seq, read, scores, matrix)
		return AlignResult{
			Alignment: alignment,
			Score:     score,
			TargetPos: refEnd - len(s.Seq) - len(seq) + targetStart,
			QueryPos:  queryStart,
			Path:      s.Path,
		}
	}

	bestResult := AlignResult{Score: math.MinInt64}
	for _, i := range n.Prev {
		currResult := LeftAlignTraversal(i.Dest, s.Seq, len(i.Dest.Seq), s.Path, extension, read, scores, matrix, pool)

		if currResult.Score > bestResult.Score {
			bestResult = currResult
		}
	}

	cigar.ReverseBytesCigar(bestResult.Alignment)
	ReversePath(bestResult.Path)

	return bestResult
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, matrix *MatrixAln) (int64, []cigar.ByteCigar, int, int) {
	alphaLen, betaLen := len(alpha), len(beta)
	route := []cigar.ByteCigar{}
	routeIdx := 0
	var i, j int

	// Reset matrix
	for i = 0; i <= alphaLen; i++ {
		matrix.m[i][0] = 0
	}
	for j = 0; j <= betaLen; j++ {
		matrix.m[0][j] = 0
	}

	for i := 1; i <= alphaLen; i++ {
		for j := 1; j <= betaLen; j++ {
			score, trace := cigar.ByteMatrixTrace(matrix.m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], matrix.m[i][j-1], matrix.m[i-1][j])

			matrix.m[i][j] = score
			matrix.trace[i][j] = trace

			if matrix.m[i][j] < 0 {
				matrix.m[i][j] = 0
			}
		}
	}

	for i, j := alphaLen, betaLen; matrix.m[i][j] > 0; {
		if len(route) == 0 || route[routeIdx].Op != matrix.trace[i][j] {
			route = append(route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[i][j]})
			routeIdx = len(route) - 1
		} else {
			route[routeIdx].RunLen++
		}

		switch matrix.trace[i][j] {
		case 'M':
			i--
			j--
		case 'I':
			j--
		case 'D':
			i--
		default:
			log.Fatalf("Error: unexpected traceback %c\n", matrix.trace[i][j])
		}
	}

	return matrix.m[alphaLen][betaLen], route, i, j
}

// 1. Reduced Memory Allocation/Deallocation
// The approach focuses on reducing unnecessary allocations and function calls.

func RightAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, extension int, read []dna.Base, scoreMatrix [][]int64, matrix *MatrixAln, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool) AlignResult {
	s := pool.Get().(*dnaPool)
	s.Seq = append(s.Seq[:0], seq[start:]...)
	s.Seq = append(s.Seq, n.Seq...)

	if len(s.Path) < len(currentPath) {
		s.Path = make([]uint32, len(currentPath))
	}
	copy(s.Path, currentPath)

	if len(seq)+len(n.Seq)-start >= extension || len(n.Next) == 0 {
		alignmentResult := RightDynamicAln(s.Seq, read, scoreMatrix, matrix, -600, dynamicScore)
		pool.Put(s)
		return AlignResult{
			Alignment: alignmentResult.Alignment,
			Score:     alignmentResult.Score,
			TargetPos: alignmentResult.TargetPos + start,
			QueryPos:  alignmentResult.QueryPos,
			Path:      s.Path,
		}
	}

	sk.rightScore = math.MinInt64
	for _, i := range n.Next {
		tempSeq := s.Seq[:len(n.Seq)]
		result := RightAlignTraversal(i.Dest, tempSeq, 0, s.Path, extension, read, scoreMatrix, matrix, sk, dynamicScore, pool)
		if result.Score > sk.rightScore {
			sk.rightScore = result.Score
			sk.rightAlignment = result.Alignment
			sk.targetEnd = result.TargetPos
			sk.queryEnd = result.QueryPos
			sk.rightPath = result.Path
		}
	}
	pool.Put(s)
	cigar.ReverseBytesCigar(sk.rightAlignment)
	return AlignResult{
		Alignment: sk.rightAlignment,
		Score:     sk.rightScore,
		TargetPos: sk.targetEnd + start,
		QueryPos:  sk.queryEnd,
		Path:      sk.rightPath,
	}
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int64, matrix *MatrixAln, gapPen int64, dynamicScore dynamicScoreKeeper) AlignResult {
	resetDynamicScore(dynamicScore)

	for i := range matrix.m {
		matrix.m[i][0] = int64(i) * gapPen
		matrix.trace[i][0] = 'D'
	}
	for j := range matrix.m[0] {
		matrix.m[0][j] = int64(j) * gapPen
		matrix.trace[0][j] = 'I'
	}

	var maxI, maxJ int
	for dynamicScore.i = 1; dynamicScore.i < len(alpha)+1; dynamicScore.i++ {
		for dynamicScore.j = 1; dynamicScore.j < len(beta)+1; dynamicScore.j++ {
			matrix.m[dynamicScore.i][dynamicScore.j], matrix.trace[dynamicScore.i][dynamicScore.j] = cigar.ByteMatrixTrace(matrix.m[dynamicScore.i-1][dynamicScore.j-1]+scores[alpha[dynamicScore.i-1]][beta[dynamicScore.j-1]], matrix.m[dynamicScore.i][dynamicScore.j-1]+gapPen, matrix.m[dynamicScore.i-1][dynamicScore.j]+gapPen)
			if matrix.m[dynamicScore.i][dynamicScore.j] > dynamicScore.currMax {
				dynamicScore.currMax = matrix.m[dynamicScore.i][dynamicScore.j]
				maxI = dynamicScore.i
				maxJ = dynamicScore.j
			}
		}
	}

	traceback(matrix, dynamicScore, maxI, maxJ)
	return AlignResult{
		Score:     matrix.m[maxI][maxJ],
		Alignment: dynamicScore.route,
		TargetPos: maxI,
		QueryPos:  maxJ,
	}
}

// Separated traceback logic for clarity
func traceback(matrix *MatrixAln, dynamicScore dynamicScoreKeeper, maxI int, maxJ int) {
	for dynamicScore.i, dynamicScore.j, dynamicScore.routeIdx = maxI, maxJ, 0; dynamicScore.i > 0 || dynamicScore.j > 0; {
		if len(dynamicScore.route) == 0 {
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[dynamicScore.i][dynamicScore.j]})
		} else if dynamicScore.route[dynamicScore.routeIdx].Op == matrix.trace[dynamicScore.i][dynamicScore.j] {
			dynamicScore.route[dynamicScore.routeIdx].RunLen += 1
		} else {
			dynamicScore.route = append(dynamicScore.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[dynamicScore.i][dynamicScore.j]})
			dynamicScore.routeIdx++
		}
		switch matrix.trace[dynamicScore.i][dynamicScore.j] {
		case 'M':
			dynamicScore.i, dynamicScore.j = dynamicScore.i-1, dynamicScore.j-1
		case 'I':
			dynamicScore.j -= 1
		case 'D':
			dynamicScore.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", matrix.trace[dynamicScore.i][dynamicScore.j])
		}
	}
}

func ReversePath(alpha []uint32) {
	n := len(alpha)
	for i := 0; i < n/2; i++ {
		alpha[i], alpha[n-i-1] = alpha[n-i-1], alpha[i]
	}
}

func leftSeed(i int) int {
	return 2*i + 1
}

func rightSeed(i int) int {
	return 2*i + 2
}

func seedsHeapify(a []SeedDev, i int) []SeedDev {
	l := leftSeed(i)
	r := rightSeed(i)
	var max int
	if l < len(a) && l >= 0 && a[l].TotalLength < a[i].TotalLength {
		max = l
	} else {
		max = i
	}
	if r < len(a) && r >= 0 && a[r].TotalLength < a[max].TotalLength {
		max = r
	}
	if max != i {
		a[i], a[max] = a[max], a[i]
		a = seedsHeapify(a, max)
	}
	return a
}

func buildSeedHeap(a []SeedDev) []SeedDev {
	for i := len(a)/2 - 1; i >= 0; i-- {
		a = seedsHeapify(a, i)
	}
	return a
}

func heapSortSeeds(a []SeedDev) {
	a = buildSeedHeap(a)
	size := len(a)
	for i := size - 1; i >= 1; i-- {
		a[0], a[i] = a[i], a[0]
		size--
		seedsHeapify(a[:size], 0)
	}
}

func quickSort(arr []*SeedDev) []*SeedDev {
	newArr := make([]*SeedDev, len(arr))
	copy(newArr, arr)
	recursiveSort(newArr, 0, len(arr)-1)
	return newArr
}

func recursiveSort(arr []*SeedDev, start, end int) {
	if (end - start) < 1 {
		return
	}

	pivot := arr[end]
	splitIndex := start

	// Iterate sub array to find values less than pivot
	// and move them to the beginning of the array
	// keeping splitIndex denoting less-value array size
	for i := start; i < end; i++ {
		if arr[i].TotalLength > pivot.TotalLength {
			if splitIndex != i {
				temp := arr[splitIndex]

				arr[splitIndex] = arr[i]
				arr[i] = temp
			}

			splitIndex++
		}
	}

	arr[end] = arr[splitIndex]
	arr[splitIndex] = pivot

	recursiveSort(arr, start, splitIndex-1)
	recursiveSort(arr, splitIndex+1, end)
}

type seedHelper struct {
	currHits                                  []uint64
	codedNodeCoord                            uint64
	seqKey                                    uint64
	keyShift                                  uint
	keyIdx, keyOffset, readOffset, nodeOffset int
	nodeIdx, nodePos                          int64
	leftMatches                               int
	rightMatches                              int
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
	if runLen < lengthOfRead {
		answer := make([]cigar.ByteCigar, 0, len(cig)+2)
		if front > 0 {
			answer = append(answer, cigar.ByteCigar{RunLen: uint16(front), Op: 'S'})
		}
		answer = append(answer, cig...)
		if front+cigar.QueryRunLen(cig) < lengthOfRead {
			answer = append(answer, cigar.ByteCigar{RunLen: uint16(lengthOfRead - front - runLen), Op: 'S'})
		}
		return answer
	} else {
		return cig
	}
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
