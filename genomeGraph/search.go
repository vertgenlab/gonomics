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

// getLeftBases retrieves bases from the left of a given position in a node's sequence and appends additional sequence.
func getLeftBases(n *Node, extension int, refEnd int, seq []dna.Base, ans []dna.Base) []dna.Base {
	seqLen := len(seq)
	basesToTake := int(math.Min(float64(seqLen+refEnd), float64(extension))) - seqLen
	if basesToTake > 0 {
		ans = append(ans, n.Seq[refEnd-basesToTake:refEnd]...)
	}
	ans = append(ans, seq...)
	return ans
}

func LeftAlignTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, read []dna.Base, settings *GraphSettings, sk scoreKeeper, memory *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	s := memory.Get().(*MatrixMemoryPool)
	defer memory.Put(s)

	s.Seq, s.Path = s.Seq[:0], s.Path[:0]
	s.Seq = getLeftBases(n, settings.Extension, refEnd, seq, s.Seq)
	s.Path = append(s.Path[:0], currentPath...)
	AddPath(s.Path, n.Id)

	if len(seq)+refEnd >= settings.Extension || len(n.Prev) == 0 {
		sk.leftScore, sk.leftAlignment, sk.targetStart, sk.queryStart = LeftDynamicAln(s.Seq, read, settings, memory)
		sk.targetStart = refEnd - len(s.Seq) - len(seq) + sk.targetStart
		sk.leftPath = s.Path
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	}

	sk.leftScore = math.MinInt64
	for _, i := range n.Prev {
		s.Route, s.CurrScore, s.TargetStart, s.QueryStart, s.Path = LeftAlignTraversal(i.Dest, s.Seq, len(i.Dest.Seq), s.Path, read, settings, sk, memory)
		if s.CurrScore > sk.leftScore {
			sk.leftScore = s.CurrScore
			sk.leftAlignment = s.Route
			sk.targetStart = s.TargetStart
			sk.queryStart = s.QueryStart
			sk.leftPath = s.Path
		}
	}
	ReversePath(sk.leftPath)
	return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, settings *GraphSettings, pool *sync.Pool) (int64, []cigar.ByteCigar, int, int) {
	rows, columns := len(alpha)+1, len(beta)+1
	matrix := pool.Get().(*MatrixMemoryPool)
	matrix.Route = matrix.Route[:0]
	matrix.CurrMax = 0

	if cap(matrix.matrix) < rows || cap(matrix.matrix[0]) < columns {
		matrix.matrix = make([][]int64, rows)
		matrix.trace = make([][]byte, rows)
		for idx := range matrix.matrix {
			matrix.matrix[idx] = make([]int64, columns)
			matrix.trace[idx] = make([]byte, columns)
		}
	}

	for matrix.i = 0; matrix.i < rows; matrix.i++ {
		matrix.matrix[matrix.i][0] = 0
	}

	for matrix.j = 0; matrix.j < columns; matrix.j++ {
		matrix.matrix[0][matrix.j] = 0
	}

	for matrix.i = 1; matrix.i < rows; matrix.i++ {
		for matrix.j = 1; matrix.j < columns; matrix.j++ {
			matrix.matrix[matrix.i][matrix.j], matrix.trace[matrix.i][matrix.j] = cigar.ByteMatrixTrace(matrix.matrix[matrix.i-1][matrix.j-1]+settings.ScoreMatrix[alpha[matrix.i-1]][beta[matrix.j-1]], matrix.matrix[matrix.i][matrix.j-1]+settings.GapPenalty, matrix.matrix[matrix.i-1][matrix.j]+settings.GapPenalty)
			if matrix.matrix[matrix.i][matrix.j] < 0 {
				matrix.matrix[matrix.i][matrix.j] = 0
			}
		}
	}

	for matrix.i, matrix.j = rows-1, columns-1; matrix.i > 0 && matrix.j > 0 && matrix.matrix[matrix.i][matrix.j] > 0; {
		switch matrix.trace[matrix.i][matrix.j] {
		case cigar.Match:
			matrix.Route = cigar.AddCigarByte(matrix.Route, cigar.ByteCigar{RunLen: 1, Op: cigar.Match})
			matrix.i--
			matrix.j--
		case cigar.Insertion:
			matrix.Route = cigar.AddCigarByte(matrix.Route, cigar.ByteCigar{RunLen: 1, Op: cigar.Insertion})
			matrix.j--
		case cigar.Deletion:
			matrix.Route = cigar.AddCigarByte(matrix.Route, cigar.ByteCigar{RunLen: 1, Op: cigar.Deletion})
			matrix.i--
		default:
			log.Fatalf("Error: unexpected traceback %c\n", matrix.trace[matrix.i][matrix.j])
		}
	}
	return matrix.matrix[rows-1][columns-1], matrix.Route, matrix.i, matrix.j
}

// getRightBases retrieves bases from the right of a given start position in a node's sequence and appends it to the provided sequence.
func getRightBases(n *Node, extension int, start int, seq []dna.Base, ans []dna.Base) []dna.Base {
	availableBases := len(n.Seq) - start
	basesToTake := int(math.Min(float64(availableBases), float64(extension)))
	ans = append(ans, seq...)
	if basesToTake > 0 {
		ans = append(ans, n.Seq[start:start+basesToTake]...)
	}
	return ans
}

func RightAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, read []dna.Base, settings *GraphSettings, sk *scoreKeeper, memory *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	s := memory.Get().(*MatrixMemoryPool)
	defer memory.Put(s)

	s.Seq, s.Path = getRightBases(n, settings.Extension, start, seq, s.Seq[:0]), append([]uint32(nil), currentPath...) // Reuse memory more efficiently
	if len(seq)+len(n.Seq)-start >= settings.Extension || len(n.Next) == 0 {
		sk.rightScore, sk.rightAlignment, sk.targetEnd, sk.queryEnd = RightDynamicAln(s.Seq, read, settings, memory)
		sk.rightPath = s.Path
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
	}

	sk.rightScore = math.MinInt64
	for _, i := range n.Next {
		alignment, score, targetEnd, queryEnd, path := RightAlignTraversal(i.Dest, s.Seq, 0, s.Path, read, settings, sk, memory)
		if score > sk.rightScore {
			sk.rightScore, sk.rightAlignment, sk.targetEnd, sk.queryEnd, sk.rightPath = score, alignment, targetEnd, queryEnd, path
		}
	}
	cigar.ReverseBytesCigar(sk.rightAlignment)
	return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, settings *GraphSettings, pool *sync.Pool) (int64, []cigar.ByteCigar, int, int) {
	rows, columns := len(alpha)+1, len(beta)+1
	matrix := pool.Get().(*MatrixMemoryPool)
	defer pool.Put(matrix)

	var maxI int
	var maxJ int

	if cap(matrix.matrix) < rows || cap(matrix.matrix[0]) < columns {
		matrix.matrix = make([][]int64, rows)
		matrix.trace = make([][]byte, rows)
		for idx := range matrix.matrix {
			matrix.matrix[idx] = make([]int64, columns)
			matrix.trace[idx] = make([]byte, columns)
		}
	}

	for matrix.i = 0; matrix.i < rows; matrix.i++ {
		for matrix.j = 0; matrix.j < columns; matrix.j++ {
			if matrix.i == 0 && matrix.j == 0 {
				matrix.matrix[matrix.i][matrix.j] = 0
			} else if matrix.i == 0 {
				matrix.matrix[matrix.i][matrix.j] = matrix.matrix[matrix.i][matrix.j-1] + settings.GapPenalty
				matrix.trace[matrix.i][matrix.j] = cigar.Insertion
			} else if matrix.j == 0 {
				matrix.matrix[matrix.i][matrix.j] = matrix.matrix[matrix.i-1][matrix.j] + settings.GapPenalty
				matrix.trace[matrix.i][matrix.j] = cigar.Deletion
			} else {
				matrix.matrix[matrix.i][matrix.j], matrix.trace[matrix.i][matrix.j] = cigar.ByteMatrixTrace(matrix.matrix[matrix.i-1][matrix.j-1]+settings.ScoreMatrix[alpha[matrix.i-1]][beta[matrix.j-1]], matrix.matrix[matrix.i][matrix.j-1]+settings.GapPenalty, matrix.matrix[matrix.i-1][matrix.j]+settings.GapPenalty)
			}
			if matrix.matrix[matrix.i][matrix.j] > matrix.CurrMax {
				matrix.CurrMax = matrix.matrix[matrix.i][matrix.j]
				maxI = matrix.i
				maxJ = matrix.j
			}
		}
	}
	for matrix.i, matrix.j, matrix.routeIdx = maxI, maxJ, 0; matrix.i > 0 || matrix.j > 0; {
		if len(matrix.Route) == 0 {
			matrix.Route = append(matrix.Route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[matrix.i][matrix.j]})
		} else if matrix.Route[matrix.routeIdx].Op == matrix.trace[matrix.i][matrix.j] {
			matrix.Route[matrix.routeIdx].RunLen += 1
		} else {
			matrix.Route = append(matrix.Route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[matrix.i][matrix.j]})
			matrix.routeIdx++
		}
		switch matrix.trace[matrix.i][matrix.j] {
		case cigar.Match:
			matrix.i, matrix.j = matrix.i-1, matrix.j-1
		case cigar.Insertion:
			matrix.j -= 1
		case cigar.Deletion:
			matrix.i -= 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", matrix.trace[matrix.i][matrix.j])
		}
	}
	return matrix.matrix[maxI][maxJ], matrix.Route, maxI, maxJ
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
