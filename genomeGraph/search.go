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

func LeftAlignTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, read []dna.Base, settings *GraphSettings, sk scoreKeeper, memory *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	cache := memory.Get().(*MatrixScore)
	defer memory.Put(cache)

	cache.Seq, cache.Path = cache.Seq[:0], cache.Path[:0]
	cache.Seq = getTargetBases(n, settings.extension, refEnd, seq, cache.Seq, left)
	cache.Path = make([]uint32, len(currentPath))
	copy(cache.Path, currentPath)
	AddPath(cache.Path, n.Id)

	currSeqLen := refEnd - len(cache.Seq) - len(seq)
	if len(seq)+refEnd >= settings.extension || len(n.Prev) == 0 {
		sk.leftScore, sk.leftAlignment, sk.targetStart, sk.queryStart = LeftDynamicAln(cache.Seq, read, settings, cache)
		sk.targetStart = currSeqLen + sk.targetStart
		sk.leftPath = cache.Path
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	} else {
		//A very negative number
		sk.leftScore = math.MinInt64
		for _, i := range n.Prev {
			cache.route, cache.currScore, cache.targetStart, cache.queryStart, cache.Path = LeftAlignTraversal(i.Dest, cache.Seq, len(i.Dest.Seq), cache.Path, read, settings, sk, memory)
			if cache.currScore > sk.leftScore {
				sk.leftScore = cache.currScore
				sk.leftAlignment = cache.route
				sk.targetStart = currSeqLen + cache.targetStart
				sk.queryStart = cache.queryStart
				sk.leftPath = cache.Path
			}
		}

		cigar.ReverseBytesCigar(sk.leftAlignment)
		ReversePath(sk.leftPath)
		return sk.leftAlignment, sk.leftScore, sk.targetStart, sk.queryStart, sk.leftPath
	}
}

func LeftDynamicAln(alpha []dna.Base, beta []dna.Base, settings *GraphSettings, matrix *MatrixScore) (int64, []cigar.ByteCigar, int, int) {
	rows, columns := len(alpha), len(beta)
	matrix.route = matrix.route[:0]
	matrix.currMax = 0

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

	for matrix.i = 1; matrix.i < rows+1; matrix.i++ {
		for matrix.j = 1; matrix.j < columns+1; matrix.j++ {
			matrix.matrix[matrix.i][matrix.j], matrix.trace[matrix.i][matrix.j] = cigar.ByteMatrixTrace(matrix.matrix[matrix.i-1][matrix.j-1]+settings.scores[alpha[matrix.i-1]][beta[matrix.j-1]], matrix.matrix[matrix.i][matrix.j-1]+settings.gapPenalty, matrix.matrix[matrix.i-1][matrix.j]+settings.gapPenalty)
			if matrix.matrix[matrix.i][matrix.j] < 0 {
				matrix.matrix[matrix.i][matrix.j] = 0
			}
		}
	}

	for matrix.i, matrix.j, matrix.routeIdx = rows, columns, 0; matrix.matrix[matrix.i][matrix.j] > 0; {
		if len(matrix.route) == 0 {
			matrix.route = append(matrix.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[matrix.i][matrix.j]})
		} else if matrix.route[matrix.routeIdx].Op == matrix.trace[matrix.i][matrix.j] {
			matrix.route[matrix.routeIdx].RunLen += 1
		} else {
			matrix.route = append(matrix.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[matrix.i][matrix.j]})
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
			log.Fatalf("Error: unexpected traceback %c\n", matrix.trace[matrix.i][matrix.j])
		}
	}
	return matrix.matrix[rows][columns], matrix.route, matrix.i, matrix.j
}

func RightAlignTraversal(n *Node, seq []dna.Base, start int, currentPath []uint32, read []dna.Base, settings *GraphSettings, sk *scoreKeeper, memory *sync.Pool) ([]cigar.ByteCigar, int64, int, int, []uint32) {
	cache := memory.Get().(*MatrixScore)
	defer memory.Put(cache)

	cache.Seq, cache.Path = cache.Seq[:0], cache.Path[:0]
	cache.Seq = getTargetBases(n, settings.extension, start, seq, cache.Seq, right)

	// Reuse the currentPath slice if possible
	if cap(cache.Path) >= len(currentPath) {
		cache.Path = cache.Path[:len(currentPath)]
	} else {
		cache.Path = make([]uint32, len(currentPath))
	}
	copy(cache.Path, currentPath)

	if len(seq)+len(n.Seq)-start >= settings.extension || len(n.Next) == 0 {
		sk.rightScore, sk.rightAlignment, sk.targetEnd, sk.queryEnd = RightDynamicAln(cache.Seq, read, settings, cache)
		sk.rightPath = cache.Path

		return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
	} else {
		sk.rightScore = math.MinInt64
		for _, i := range n.Next {
			cache.route, cache.currScore, cache.targetEnd, cache.queryEnd, cache.Path = RightAlignTraversal(i.Dest, cache.Seq, 0, cache.Path, read, settings, sk, memory)
			if cache.currScore > sk.rightScore {
				sk.rightScore = cache.currScore
				sk.rightAlignment = cache.route
				sk.targetEnd = cache.targetEnd
				sk.queryEnd = cache.queryEnd
				sk.rightPath = cache.Path
			}
		}
		cigar.ReverseBytesCigar(sk.rightAlignment)
		return sk.rightAlignment, sk.rightScore, sk.targetEnd + start, sk.queryEnd, sk.rightPath
	}
}

func RightDynamicAln(alpha []dna.Base, beta []dna.Base, settings *GraphSettings, matrix *MatrixScore) (int64, []cigar.ByteCigar, int, int) {
	rows, columns := len(alpha)+1, len(beta)+1
	matrix.route = matrix.route[:0]
	matrix.currMax = 0

	if cap(matrix.matrix) < rows || cap(matrix.matrix[0]) < columns {
		matrix.matrix = make([][]int64, rows)
		matrix.trace = make([][]byte, rows)
		for idx := range matrix.matrix {
			matrix.matrix[idx] = make([]int64, columns)
			matrix.trace[idx] = make([]byte, columns)
		}
	}
	var maxI, maxJ int = 0, 0

	for matrix.j = 0; matrix.j < columns; matrix.j++ {
		matrix.matrix[0][matrix.j] = int64(matrix.j) * settings.gapPenalty
		matrix.trace[0][matrix.j] = cigar.Insertion
	}
	for matrix.i = 0; matrix.i < rows; matrix.i++ {
		matrix.matrix[matrix.i][0] = int64(matrix.i) * settings.gapPenalty
		matrix.trace[matrix.i][0] = cigar.Deletion
	}

	for matrix.i = 1; matrix.i < rows; matrix.i++ {
		for matrix.j = 1; matrix.j < columns; matrix.j++ {
			matrix.matrix[matrix.i][matrix.j], matrix.trace[matrix.i][matrix.j] = cigar.ByteMatrixTrace(matrix.matrix[matrix.i-1][matrix.j-1]+settings.scores[alpha[matrix.i-1]][beta[matrix.j-1]], matrix.matrix[matrix.i][matrix.j-1]+settings.gapPenalty, matrix.matrix[matrix.i-1][matrix.j]+settings.gapPenalty)
			if matrix.matrix[matrix.i][matrix.j] > matrix.currMax {
				matrix.currMax = matrix.matrix[matrix.i][matrix.j]
				maxI = matrix.i
				maxJ = matrix.j
			}
		}
	}

	for matrix.i, matrix.j, matrix.routeIdx = maxI, maxJ, 0; matrix.i > 0 || matrix.j > 0; {
		if len(matrix.route) == 0 {
			matrix.route = append(matrix.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[matrix.i][matrix.j]})
		} else if matrix.route[matrix.routeIdx].Op == matrix.trace[matrix.i][matrix.j] {
			matrix.route[matrix.routeIdx].RunLen += 1
		} else {
			matrix.route = append(matrix.route, cigar.ByteCigar{RunLen: 1, Op: matrix.trace[matrix.i][matrix.j]})
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
	return matrix.matrix[maxI][maxJ], matrix.route, maxI, maxJ
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

func seedsHeapify(a []SeedDev, n, i int) {
	largest := i
	l := leftSeed(i)
	r := rightSeed(i)

	if l < n && a[l].TotalLength < a[largest].TotalLength {
		largest = l
	}
	if r < n && a[r].TotalLength < a[largest].TotalLength {
		largest = r
	}
	if largest != i {
		a[i], a[largest] = a[largest], a[i]
		seedsHeapify(a, n, largest)
	}
}

func buildSeedHeap(a []SeedDev) {
	n := len(a)
	for i := n/2 - 1; i >= 0; i-- {
		seedsHeapify(a, n, i)
	}
}

func heapSortSeeds(a []SeedDev) {
	buildSeedHeap(a)
	n := len(a)
	for i := n - 1; i > 0; i-- {
		a[0], a[i] = a[i], a[0]
		seedsHeapify(a, i, 0)
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
func seedMapMemPool(seedHash map[uint64][]uint64, nodes []Node, read *fastq.FastqBig, seedLen int, finalSeeds []SeedDev, tempSeeds []SeedDev, seedBuildHelper *seedHelper) []SeedDev {
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
