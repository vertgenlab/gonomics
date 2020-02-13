package simpleGraph

import (
	//"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	//"sync"
)

func AlignTraversalFwd(n *Node, seq []dna.Base, start int, currentPath []uint32, extention int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, []uint32) {
	currentPath = AddPath(n.Id, currentPath)
	var bestQueryEnd, queryEnd int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var path, bestPath []uint32

	if len(seq) >= extention {
		log.Fatalf("Error: the length of DNA sequence in previous nodes should not be enough to satisfy the desired extenion.\n")
	}
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, extention)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)

	var base int
	for base = 0; base < len(seq); base++ {
		s[base] = seq[base]
	}
	for base = 0; base < len(n.Seq[start:start+basesToTake]); base++ {
		s[len(seq)+base] = n.Seq[start+base]
	}
	if availableBases >= extention || len(n.Next) == 0 {
		score, alignment, _, _, _, queryEnd = RightLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, queryEnd, currentPath
	} else {

		for _, i := range n.Next {
			bestScore = -1
			alignment, score, queryEnd, path = AlignTraversalFwd(i.Dest, s, 0, currentPath, extention, read, m, trace)
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

func AlignReverseGraphTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, extention int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, int, []uint32) {
	currentPath = AddPath(n.Id, currentPath)
	var bestQueryStart, queryStart, refStart, bestRefStart int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var path, bestPath []uint32

	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, extention)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	var base int
	for base = 0; base < targetLength; base++ {
		if base < basesToTake {
			s[base] = n.Seq[refEnd-basesToTake+base]
		} else {
			s[base] = seq[base]
		}
	}
	if availableBases >= extention || len(n.Next) == 0 {
		score, alignment, refStart, _, queryStart, _ = LeftLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, refEnd - basesToTake + refStart, queryStart, currentPath
	} else {
		bestScore = -1
		//tmp := make([]uint32, len(currentPath))
		//copy(tmp, currentPath)
		for _, i := range n.Prev {
			alignment, score, refStart, queryStart, path = AlignReverseGraphTraversal(i.Dest, s, len(i.Dest.Seq), currentPath, extention, read, m, trace)
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
