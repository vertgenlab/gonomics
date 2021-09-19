package genomeGraph

import (
	"bytes"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
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
	m     [][]int
	trace [][]byte
}

type dnaPool struct {
	Seq         []dna.Base
	Path        []uint32
	queryStart  int
	queryEnd    int
	targetStart int
	targetEnd   int
	currScore   int
}

type dynamicScoreKeeper struct {
	i        int
	j        int
	routeIdx int
	currMax  int
	route    []cigar.ByteCigar
}

type scoreKeeper struct {
	targetStart  int
	targetEnd    int
	queryStart   int
	queryEnd     int
	extension    int
	currScore    int
	seedScore    int
	perfectScore int
	leftScore    int
	rightScore   int
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

func MatrixSetup(size int) ([][]int, [][]byte) {
	m := make([][]int, size)
	trace := make([][]byte, size)
	for idx := range m {
		m[idx] = make([]int, size)
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

func getLeftTargetBases(n *Node, extension int, refEnd int, seq []dna.Base, ans []dna.Base) []dna.Base {
	//var availableBases int = len(seq) + refEnd
	//var targetLength int = common.Min(availableBases, extension)
	//var basesToTake int = targetLength - len(seq)
	return append(append(ans, n.Seq[refEnd-numbers.Min(len(seq)+refEnd, extension)-len(seq):refEnd]...), seq...)
}

func getRightBases(n *Node, extension int, start int, seq []dna.Base, ans []dna.Base) []dna.Base {
	//var availableBases int = len(seq) + len(n.Seq) - start
	//var targetLength int = common.Min(availableBases, extension)
	//var basesToTake int = targetLength - len(seq)
	return append(append(ans, seq...), n.Seq[start:start+numbers.Min(len(seq)+len(n.Seq)-start, extension)-len(seq)]...)
}

/*
func leftBasesFromTwoBit(n *Node, extension int, refEnd int, seq []dna.Base, ans []dna.Base) []dna.Base {
	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, extension)
	var basesToTake int = targetLength - len(seq)
	return append(append(ans, seq...), dnaTwoBit.GetFrag(n.SeqTwoBit, refEnd-basesToTake, refEnd)...)
}*/

/*
func rightBasesFromTwoBit(n *Node, extension int, start int, seq []dna.Base, ans []dna.Base) []dna.Base {
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, extension)
	var basesToTake int = targetLength - len(seq)

	return append(append(ans, seq...), dnaTwoBit.GetFrag(n.SeqTwoBit, start, start+basesToTake)...)
}*/

func LeftAlignTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, extension int, read []dna.Base, scores [][]int, matrix *MatrixAln, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool) ([]cigar.ByteCigar, int, int, int, []uint32) {
	//if len(seq) >= extension {
	//	log.Fatalf("Error: left traversal, the length=%d of DNA sequence in previous nodes should not be enough to satisfy the desired extenion=%d.\n", len(seq), extension)
	//}
	s := pool.Get().(*dnaPool)
	s.Seq, s.Path = s.Seq[:0], s.Path[:0]
	s.Seq = getLeftTargetBases(n, extension, refEnd, seq, s.Seq)
	s.Path = make([]uint32, len(currentPath))
	copy(s.Path, currentPath)
	AddPath(s.Path, n.Id)
	if len(seq)+refEnd >= extension || len(n.Prev) == 0 {
		sk.leftScore, sk.leftAlignment, sk.targetStart, sk.queryStart = LeftDynamicAln(s.Seq, read, scores, matrix, -600, dynamicScore)
		sk.targetStart = refEnd - len(s.Seq) - len(seq) + sk.targetStart
		sk.leftPath = s.Path
		pool.Put(s)
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	} else {
		//A very negative number
		sk.leftScore = -9223372036854775808 //TODO: fix math.MinInt
		for _, i := range n.Prev {
			dynamicScore.route, s.currScore, s.targetStart, s.queryStart, s.Path = LeftAlignTraversal(i.Dest, s.Seq, len(i.Dest.Seq), s.Path, extension, read, scores, matrix, sk, dynamicScore, pool)
			if s.currScore > sk.leftScore {
				sk.leftScore = s.currScore
				sk.leftAlignment = dynamicScore.route
				sk.targetStart = refEnd - len(s.Seq) - len(seq) + s.targetStart
				sk.queryStart = s.queryStart
				sk.leftPath = s.Path
			}
		}
		pool.Put(s)
		cigar.ReverseBytesCigar(sk.leftAlignment)
		ReversePath(sk.leftPath)
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	}
}

func RightAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, extension int, read []dna.Base, scoreMatrix [][]int, matrix *MatrixAln, sk scoreKeeper, dynamicScore dynamicScoreKeeper, pool *sync.Pool) ([]cigar.ByteCigar, int, int, int, []uint32) {
	//if len(seq) >= extension {
	//	log.Fatalf("Error: right traversal, the length=%d of DNA sequence in previous nodes should not be enough to satisfy the desired extenion=%d.\n", len(seq), extension)
	//}
	s := pool.Get().(*dnaPool)
	s.Seq, s.Path = s.Seq[:0], s.Path[:0]
	s.Seq = getRightBases(n, extension, start, seq, s.Seq)
	s.Path = make([]uint32, len(currentPath))
	copy(s.Path, currentPath)
	if len(seq)+len(n.Seq)-start >= extension || len(n.Next) == 0 {
		sk.rightScore, sk.rightAlignment, sk.targetEnd, sk.queryEnd = RightDynamicAln(s.Seq, read, scoreMatrix, matrix, -600, dynamicScore)
		sk.rightPath = s.Path
		pool.Put(s)
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
	} else {
		sk.rightScore = -9223372036854775808 //TODO: fix math.MinInt64
		for _, i := range n.Next {
			dynamicScore.route, s.currScore, s.targetEnd, s.queryEnd, s.Path = RightAlignTraversal(i.Dest, s.Seq, 0, s.Path, extension, read, scoreMatrix, matrix, sk, dynamicScore, pool)
			if s.currScore > sk.rightScore {
				sk.rightScore = s.currScore
				sk.rightAlignment = dynamicScore.route
				sk.targetEnd = s.targetEnd
				sk.queryEnd = s.queryEnd
				sk.rightPath = s.Path
			}
		}
		pool.Put(s)
		cigar.ReverseBytesCigar(sk.rightAlignment)
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
	}
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int, matrix *MatrixAln, gapPen int, dynamicScore dynamicScoreKeeper) (int, []cigar.ByteCigar, int, int) {
	resetDynamicScore(dynamicScore)
	for dynamicScore.i = 0; dynamicScore.i < len(alpha)+1; dynamicScore.i++ {
		matrix.m[dynamicScore.i][0] = 0
	}

	for dynamicScore.j = 0; dynamicScore.j < len(beta)+1; dynamicScore.j++ {
		matrix.m[0][dynamicScore.j] = 0
	}

	for dynamicScore.i = 1; dynamicScore.i < len(alpha)+1; dynamicScore.i++ {
		for dynamicScore.j = 1; dynamicScore.j < len(beta)+1; dynamicScore.j++ {
			matrix.m[dynamicScore.i][dynamicScore.j], matrix.trace[dynamicScore.i][dynamicScore.j] = cigar.ByteMatrixTrace(matrix.m[dynamicScore.i-1][dynamicScore.j-1]+scores[alpha[dynamicScore.i-1]][beta[dynamicScore.j-1]], matrix.m[dynamicScore.i][dynamicScore.j-1]+gapPen, matrix.m[dynamicScore.i-1][dynamicScore.j]+gapPen)
			if matrix.m[dynamicScore.i][dynamicScore.j] < 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = 0
			}
		}
	}

	for dynamicScore.i, dynamicScore.j, dynamicScore.routeIdx = len(alpha), len(beta), 0; matrix.m[dynamicScore.i][dynamicScore.j] > 0; {
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
			log.Fatalf("Error: unexpected traceback %c\n", matrix.trace[dynamicScore.i][dynamicScore.j])
		}
	}
	return matrix.m[len(alpha)][len(beta)], dynamicScore.route, dynamicScore.i, dynamicScore.j
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, scores [][]int, matrix *MatrixAln, gapPen int, dynamicScore dynamicScoreKeeper) (int, []cigar.ByteCigar, int, int) {
	resetDynamicScore(dynamicScore)
	var maxI int
	var maxJ int
	for dynamicScore.i = 0; dynamicScore.i < len(alpha)+1; dynamicScore.i++ {
		for dynamicScore.j = 0; dynamicScore.j < len(beta)+1; dynamicScore.j++ {
			if dynamicScore.i == 0 && dynamicScore.j == 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = 0
			} else if dynamicScore.i == 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = matrix.m[dynamicScore.i][dynamicScore.j-1] + gapPen
				matrix.trace[dynamicScore.i][dynamicScore.j] = 'I'
			} else if dynamicScore.j == 0 {
				matrix.m[dynamicScore.i][dynamicScore.j] = matrix.m[dynamicScore.i-1][dynamicScore.j] + gapPen
				matrix.trace[dynamicScore.i][dynamicScore.j] = 'D'
			} else {
				matrix.m[dynamicScore.i][dynamicScore.j], matrix.trace[dynamicScore.i][dynamicScore.j] = cigar.ByteMatrixTrace(matrix.m[dynamicScore.i-1][dynamicScore.j-1]+scores[alpha[dynamicScore.i-1]][beta[dynamicScore.j-1]], matrix.m[dynamicScore.i][dynamicScore.j-1]+gapPen, matrix.m[dynamicScore.i-1][dynamicScore.j]+gapPen)
			}
			if matrix.m[dynamicScore.i][dynamicScore.j] > dynamicScore.currMax {
				dynamicScore.currMax = matrix.m[dynamicScore.i][dynamicScore.j]
				maxI = dynamicScore.i
				maxJ = dynamicScore.j
			}
		}
	}
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
	return matrix.m[maxI][maxJ], dynamicScore.route, maxI, maxJ
}

func ReversePath(alpha []uint32) {
	var i, off int
	for i = len(alpha)/2 - 1; i >= 0; i-- {
		off = len(alpha) - 1 - i
		alpha[i], alpha[off] = alpha[off], alpha[i]
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

	for i, v := range arr {
		newArr[i] = v
	}
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
	//   and move them to the beginning of the array
	//   keeping splitIndex denoting less-value array size
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
	nodeIdx, nodePos                          int
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

//seedBuildHelper.nodeIdx, seedBuildHelper.nodePos int = 0, 0
func seedMapMemPool(seedHash map[uint64][]uint64, nodes []Node, read *fastq.FastqBig, seedLen int, perfectScore int, scoreMatrix [][]int, finalSeeds []SeedDev, tempSeeds []SeedDev, seedBuildHelper *seedHelper) []SeedDev {
	const basesPerInt int = 32
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
		common.ExitIfError(err)
		err = buf.WriteByte('\n')
		common.ExitIfError(err)
		_, err = buf.WriteString(giraf.ToString(&gp.Rev))
		common.ExitIfError(err)
		err = buf.WriteByte('\n')
		common.ExitIfError(err)
		io.Copy(file, buf)
		simplePool.Put(buf)
	}
	file.Close()
	wg.Done()
}
