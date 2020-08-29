package simpleGraph

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
)

func leftSeed(i int) int {
	return 2*i + 1
}

func rightSeed(i int) int {
	return 2*i + 2
}

func seedsHeapify(a []*SeedDev, i int) []*SeedDev {
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

func buildSeedHeap(a []*SeedDev) []*SeedDev {
	for i := len(a)/2 - 1; i >= 0; i-- {
		a = seedsHeapify(a, i)
	}
	return a
}

func heapSortSeeds(a []*SeedDev) {
	a = buildSeedHeap(a)
	size := len(a)
	for i := size - 1; i >= 1; i-- {
		a[0], a[i] = a[i], a[0]
		size--
		seedsHeapify(a[:size], 0)
	}
}

func extendToTheRightDev(node *Node, read *fastq.FastqBig, readStart int, nodeStart int, posStrand bool, answer []*SeedDev) []*SeedDev {
	const basesPerInt int = 32
	answer = answer[:0]
	var nodeOffset int = nodeStart % basesPerInt
	var readOffset int = 31 - ((readStart - nodeOffset + 31) % 32)
	var rightMatches int = 0
	var currNode SeedDev
	var nextParts []*SeedDev
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
			//nextParts := extendSeeds.Get().([]*SeedDev)
			nextParts = extendToTheRightDev(node.Next[i].Dest, read, readStart+rightMatches, 0, posStrand, nextParts)
			// if we aligned into the next node, make a seed for this node and point it to the next one

			for j = 0; j < len(nextParts); j++ {
				currNode = SeedDev{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(rightMatches), PosStrand: posStrand, TotalLength: uint32(rightMatches) + nextParts[j].TotalLength, NextPart: nextParts[j], Next: nil}
				answer = append(answer, &currNode)
			}

			//nextParts = nextParts[:0]
			//extendSeeds.Put(nextParts)
		}
	}
	// if the alignment did not go to another node, return the match for this node
	if len(answer) == 0 {
		currNode = SeedDev{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(rightMatches), PosStrand: posStrand, TotalLength: uint32(rightMatches), NextPart: nil, Next: nil}
		answer = []*SeedDev{&currNode}
	}

	return answer
}

func extendToTheLeftDev(node *Node, read *fastq.FastqBig, currPart *SeedDev) []*SeedDev {
	var answer, prevParts []*SeedDev
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
		return []*SeedDev{currPart}
	} else {
		return answer
	}
}

func extendToTheLeftHelperDev(node *Node, read *fastq.FastqBig, nextPart *SeedDev) []*SeedDev {
	const basesPerInt int = 32
	var nodePos int = node.SeqTwoBit.Len - 1
	var readPos int = int(nextPart.QueryStart) - 1
	var nodeOffset int = nodePos % basesPerInt
	var readOffset int = 31 - ((readPos - nodeOffset + 31) % 32)
	var leftMatches int = 0
	var currPart *SeedDev = nil
	var prevParts, answer []*SeedDev
	var i int
	var readBase dna.Base

	if nextPart.PosStrand {
		leftMatches = common.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.Rainbow[readOffset], readPos+readOffset))
	} else {
		leftMatches = common.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.RainbowRc[readOffset], readPos+readOffset))
	}

	if leftMatches == 0 {
		log.Fatal("Error: should not have zero matches to the left\n")
	}

	currPart = &SeedDev{TargetId: node.Id, TargetStart: uint32(nodePos - (leftMatches - 1)), QueryStart: uint32(readPos - (leftMatches - 1)), Length: uint32(leftMatches), PosStrand: nextPart.PosStrand, TotalLength: uint32(leftMatches) + nextPart.TotalLength, NextPart: nextPart, Next: nil}
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
		answer = []*SeedDev{currPart}
	}
	return answer
}

func seedMapMemPool(seedHash map[uint64][]uint64, nodes []*Node, read *fastq.FastqBig, seedLen int, perfectScore int64, scoreMatrix [][]int64, finalSeeds []*SeedDev, tempSeeds []*SeedDev) []*SeedDev {
	const basesPerInt int64 = 32
	var currHits []uint64
	var codedNodeCoord uint64
	var leftMatches int = 0
	var keyIdx, keyOffset, readOffset, nodeOffset int = 0, 0, 0, 0
	var nodeIdx, nodePos int64 = 0, 0
	//var poolHead *SeedDev = *memoryPool
	var seqKey uint64
	var keyShift uint = 64 - (uint(seedLen) * 2)
	//var tempSeeds []*SeedDev
	//tempSeeds := seedPool.Get().([]*SeedDev)
	var tempSeed *SeedDev
	for readStart := 0; readStart < len(read.Seq)-seedLen+1; readStart++ {
		keyIdx = (readStart + 31) / 32
		keyOffset = 31 - ((readStart + 31) % 32)

		// do fwd strand
		seqKey = read.Rainbow[keyOffset].Seq[keyIdx] >> keyShift
		currHits = seedHash[seqKey]

		for _, codedNodeCoord = range currHits {
			nodeIdx, nodePos = numberToChromAndPos(codedNodeCoord)
			nodeOffset = int(nodePos % basesPerInt)
			readOffset = 31 - ((readStart - nodeOffset + 31) % 32)

			leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), &read.Rainbow[readOffset], readStart+readOffset))
			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.Rainbow[readOffset], readStart+readOffset)
			tempSeeds = extendToTheRightDev(nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), true, tempSeeds)
			//log.Printf("After extendToTheRightDev fwd:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				finalSeeds = append(finalSeeds, extendToTheLeftDev(nodes[nodeIdx], read, tempSeed)...)
			}
		}
		// do rev strand
		seqKey = read.RainbowRc[keyOffset].Seq[keyIdx] >> keyShift
		currHits = seedHash[seqKey]
		/*if len(currHits) > 0 {
			log.Printf(" %d hits\n", len(currHits))
		}*/
		for _, codedNodeCoord = range currHits {
			nodeIdx, nodePos = numberToChromAndPos(codedNodeCoord)
			nodeOffset = int(nodePos % basesPerInt)
			readOffset = 31 - ((readStart - nodeOffset + 31) % 32)

			leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), &read.RainbowRc[readOffset], readStart+readOffset))
			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.RainbowRc[readOffset], readStart+readOffset)
			tempSeeds = extendToTheRightDev(nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), false, tempSeeds)
			//log.Printf("After extendToTheRightDev rev:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				//log.Printf("tempSeed.QueryStart = %d\n", tempSeed.QueryStart)
				finalSeeds = append(finalSeeds, extendToTheLeftDev(nodes[nodeIdx], read, tempSeed)...)
			}

		}
	}
	return finalSeeds
}
