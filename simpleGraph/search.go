package simpleGraph

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

func AlignTraversalFwd(n *Node, seq []dna.Base, start int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, []uint32) {
	currentPath = append(currentPath, n.Id)
	var bestQueryEnd, queryEnd int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var path, bestPath []uint32

	if len(seq) >= ext {
		log.Fatalf("Error: the length of DNA sequence in previous nodes should not be enough to satisfy the desired extenion.\n")
	}
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	//log.Printf("len(seq)=%d, len(n.Seq)=%d, start=%d, targetLength=%d, basesToTake=%d\n", len(seq), len(n.Seq), start, targetLength, basesToTake)
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):targetLength], n.Seq[start:start+basesToTake])

	if availableBases >= ext || len(n.Next) == 0 {
		score, alignment, _, _, _, queryEnd = RightLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, queryEnd, currentPath
	} else {
		bestScore = -1
		for _, i := range n.Next {
			tmpPath := make([]uint32, len(currentPath)) //TODO: should not need to copy path N times, but N-1 times
			copy(tmpPath, currentPath)
			alignment, score, queryEnd, path = AlignTraversalFwd(i.Dest, s, 0, tmpPath, ext, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestQueryEnd = queryEnd
				bestPath = path
			}
		}
		return bestAlignment, bestScore, bestQueryEnd, bestPath
	}
}

func AlignReverseGraphTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, int, []uint32) {
	currentPath = append([]uint32{n.Id}, currentPath...)
	var bestQueryStart, queryStart, refStart, bestRefStart int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var path, bestPath []uint32

	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	copy(s[0:basesToTake], n.Seq[refEnd-basesToTake:refEnd])
	copy(s[basesToTake:targetLength], seq)

	//log.Printf("left(reverse) alignment: seq1=%s, seq2=%s\n", dna.BasesToString(s), dna.BasesToString(read))
	if availableBases >= ext || len(n.Next) == 0 {
		//log.Printf("at leaf, about to align, path is:%v\n", currentPath)
		score, alignment, refStart, _, queryStart, _ = LeftLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, refEnd - basesToTake + refStart, queryStart, currentPath
	} else {
		bestScore = -1
		tmp := make([]uint32, len(currentPath))
		copy(tmp, currentPath)
		for _, i := range n.Prev {

			alignment, score, refStart, queryStart, path = AlignReverseGraphTraversal(i.Dest, s, len(i.Dest.Seq), currentPath, ext, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestRefStart = refStart
				bestQueryStart = queryStart
				bestPath = path
			}
		}
		return bestAlignment, bestScore, bestRefStart, bestQueryStart, bestPath
	}
}

func getSeqTraversal(curr *Node, seq []dna.Base, start int, extension int) [][]dna.Base {
	var answer [][]dna.Base
	if len(seq) >= extension {
		log.Fatalf("Error: the length of DNA sequence in previous nodes should not be enough to satisfy the desired extension.\n")
	}
	var availableBases int = len(curr.Seq) - start + len(seq)
	var targetLength int = common.Min(availableBases, extension)
	var basesToTake int = targetLength - len(seq)
	//log.Printf("len(seq)=%d, len(n.Seq)=%d, start=%d, targetLength=%d, basesToTake=%d\n", len(seq), len(curr.Seq), start, targetLength, basesToTake)
	var s []dna.Base = make([]dna.Base, targetLength)
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):targetLength], curr.Seq[start:start+basesToTake])
	if availableBases >= extension || len(curr.Next) == 0 {
		if dna.CountBaseInterval(s, dna.N, 0, len(s)) == 0 {
			answer = append(answer, s)
		}

		//score, alignment, _, _, _, queryEnd = RightLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return answer
	} else {

		for _, i := range curr.Next {
			answer = append(answer, getSeqTraversal(i.Dest, s, 0, extension)...)
		}
		return answer
	}
}

func devIndexGraph(genome []*Node, seedLen int, seedStep int) map[uint64][]*SeedBed {
	if seedLen < 2 || seedLen > 32 {
		log.Fatalf("Error: seed length needs to be greater than 1 and less than 33.  Got: %d\n", seedLen)
	}
	answer := make(map[uint64][]*SeedBed)
	var traversalSeq [][]dna.Base
	var seqCode uint64
	var nodeIdx, pos, i int
	for nodeIdx = 0; nodeIdx < len(genome); nodeIdx++ {
		for pos = 0; pos < len(genome[nodeIdx].Seq)-seedLen+1; pos += seedStep {
			traversalSeq = getSeqTraversal(genome[nodeIdx], []dna.Base{}, pos, seedLen)
			for i = 0; i < len(traversalSeq); i++ {
				seqCode = dnaToNumber(traversalSeq[i], 0, len(traversalSeq[i]))

				curr := SeedBed{Id: genome[nodeIdx].Id, Start: uint32(pos), End: uint32(pos + common.Min(len(traversalSeq[i]), seedLen)), Next: nil}
				answer[seqCode] = append(answer[seqCode], &curr)
			}
		}
	}
	return answer
}

func PointToBases(currSeq []dna.Base, incomingSeq []dna.Base, currStart int, incomingStart int, currEnd int, incomingEnd int, currFront bool) {
	newSize := currEnd + incomingEnd - currStart - incomingStart
	var i int
	if len(currSeq) < newSize {
		currSeq = append(currSeq, make([]dna.Base, newSize-len(currSeq))...)
	}
	if currFront {
		for i = currEnd; i < len(incomingSeq[incomingStart:incomingEnd]); i++ {
			currSeq[i] = incomingSeq[incomingStart+i]
		}
	}
	if !currFront {
		for i = 0; i < len(currSeq); i++ {
			currSeq[incomingEnd-incomingStart-1] = currSeq[i]
		}
		for i = 0; i < len(incomingSeq[incomingStart:incomingEnd]); i++ {
			currSeq[i] = incomingSeq[incomingStart+i]
		}
	}
}
