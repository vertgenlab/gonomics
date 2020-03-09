package simpleGraph

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

func AlignTraversalFwd(rightNode *Node, seq []dna.Base, start int, currentPath []uint32, extention int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, []uint32) {

	var bestQueryEnd, queryEnd int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var bestPath []uint32

	if len(seq) >= extention {
		log.Fatalf("Error: the length of DNA sequence in previous nodes should not be enough to satisfy the desired extension.\n")
	}
	var availableBases int = len(rightNode.Seq) - start + len(seq)

	var targetLength int = common.Min(availableBases, extention)
	var basesToTake int = targetLength - len(seq)
	//log.Printf("len(seq)=%d, len(n.Seq)=%d, start=%d, targetLength=%d, basesToTake=%d\n", len(seq), len(rightNode.Seq), start, targetLength, basesToTake)
	var s []dna.Base = make([]dna.Base, targetLength)

	copy(s[0:len(seq)], seq)
	copy(s[len(seq):targetLength], rightNode.Seq[start:start+basesToTake])

	/*
		var base int
		for base = 0; base < len(seq); base++ {
			s[base] = seq[base]
		}
		for base = 0; base < len(n.Seq[start:start+basesToTake]); base++ {
			s[len(seq)+base] = n.Seq[start+base]
		}*/
	if availableBases >= extention || len(rightNode.Next) == 0 {
		currentPath = AddPath(rightNode.Id, currentPath)
		score, alignment, _, _, _, queryEnd = RightLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, queryEnd, currentPath
	} else {

		for _, i := range rightNode.Next {
			bestScore = -1
			alignment, score, queryEnd, currentPath = AlignTraversalFwd(i.Dest, s, 0, currentPath, extention, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestQueryEnd = queryEnd
				bestPath = currentPath
			}
		}
		return bestAlignment, bestScore, bestQueryEnd, bestPath
	}
}

func AlignReverseGraphTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, extention int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, int, []uint32) {

	var bestQueryStart, queryStart, refStart, bestRefStart int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var bestPath []uint32

	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, extention)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	copy(s[0:basesToTake], n.Seq[refEnd-basesToTake:refEnd])
	copy(s[basesToTake:targetLength], seq)
	/*var base int
	for base = 0; base < targetLength; base++ {
		if base < basesToTake {
			s[base] = n.Seq[refEnd-basesToTake+base]
		} else {
			s[base] = seq[base]
		}
	}*/
	if availableBases >= extention || len(n.Next) == 0 {
		currentPath = AddPath(n.Id, currentPath)
		score, alignment, refStart, _, queryStart, _ = LeftLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, refEnd - basesToTake + refStart, queryStart, currentPath
	} else {
		bestScore = -1
		//tmp := make([]uint32, len(currentPath))
		//copy(tmp, currentPath)
		for _, i := range n.Prev {
			alignment, score, refStart, queryStart, currentPath = AlignReverseGraphTraversal(i.Dest, s, len(i.Dest.Seq), currentPath, extention, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestRefStart = refStart
				bestQueryStart = queryStart
				bestPath = currentPath
			}
		}
		reversePath(bestPath)
		return bestAlignment, bestScore, bestRefStart, bestQueryStart, bestPath
	}
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
