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

// AlignGraphTraversal performs graph traversal for alignment either to the left or right, calculating optimal alignment using dynamic programming.
func AlignGraphTraversal(n *Node, seq []dna.Base, position int, currentPath []uint32, extension int, read []dna.Base, scores [][]int64, matrix *MatrixAln, dynamicScore *dynamicScoreKeeper, pool *sync.Pool, direction byte) ([]cigar.ByteCigar, int64, int, int, []uint32) {
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
		s.Path = append([]uint32(nil), s.Path...)
		pool.Put(s) // Return the object to the pool after copying the path to avoid data races.
		return alignment, score, targetStart, queryStart, s.Path
	}

	// Recursive case: traverse next or previous nodes based on direction.
	score = math.MinInt64
	nodes := n.Next
	if direction == leftTraversal {
		nodes = n.Prev
	}
	for _, edge := range nodes {
		route, currScore, end, start, path := AlignGraphTraversal(edge.Dest, s.Seq, 0, s.Path, extension, read, scores, matrix, dynamicScore, pool, direction)
		if currScore > score {
			// Update alignment details if current score is better.
			score, alignment, targetStart, queryStart = currScore, route, end, start
			s.Path = append([]uint32(nil), path...)
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
	var currOp byte
	alignment := make([]cigar.ByteCigar, 0)
	i, j = dynamicScore.i, dynamicScore.j
	for i > 0 && j > 0 {
		switch matrix.trace[i][j] {
		case cigar.Match:
			if currOp == cigar.Match {
				alignment[len(alignment)-1].RunLen++
			} else {
				alignment = append(alignment, cigar.ByteCigar{RunLen: 1, Op: cigar.Match})
				currOp = cigar.Match
			}
			i--
			j--
		case cigar.Insertion:
			if currOp == cigar.Insertion {
				alignment[len(alignment)-1].RunLen++
			} else {
				alignment = append(alignment, cigar.ByteCigar{RunLen: 1, Op: cigar.Insertion})
				currOp = cigar.Insertion
			}
			j--
		case cigar.Deletion:
			if currOp == cigar.Deletion {
				alignment[len(alignment)-1].RunLen++
			} else {
				alignment = append(alignment, cigar.ByteCigar{RunLen: 1, Op: cigar.Deletion})
				currOp = cigar.Deletion
			}
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
