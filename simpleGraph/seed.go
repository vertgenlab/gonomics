package simpleGraph

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
)

type Seed struct {
	TargetId    uint32
	TargetStart int32
	QueryStart  int32
	Length      int32
	PosStrand   bool
	Next        *Seed
	Prev        *Seed
}

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
	nodeIdx, nodePos                          int64
	leftMatches                               int
	rightMatches                              int
	tempSeed                                  *SeedDev
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
			nextParts = extendToTheRightDev(node.Next[i].Dest, read, readStart+rightMatches, 0, posStrand, nextParts)
			// if we aligned into the next node, make a seed for this node and point it to the next one
			for j = 0; j < len(nextParts); j++ {
				currNode = SeedDev{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(rightMatches), PosStrand: posStrand, TotalLength: uint32(rightMatches) + nextParts[j].TotalLength, NextPart: nextParts[j], Next: nil}
				answer = append(answer, &currNode)
			}
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

func newSeedBuilder() *seedHelper {
	var tmp SeedDev = SeedDev{}
	return &seedHelper{
		currHits: make([]uint64, 0, 20),
		tempSeed: &tmp,
	}
}

func restartSeedHelper(helper *seedHelper) {
	helper.currHits = helper.currHits[:0]
	helper.keyIdx, helper.keyOffset, helper.readOffset, helper.nodeOffset = 0, 0, 0, 0
	helper.nodeIdx, helper.nodePos = 0, 0
	helper.seqKey, helper.codedNodeCoord = 0, 0
	helper.leftMatches = 0
}

//seedBuildHelper.nodeIdx, seedBuildHelper.nodePos int64 = 0, 0
func seedMapMemPool(seedHash map[uint64][]uint64, nodes []*Node, read *fastq.FastqBig, seedLen int, perfectScore int64, scoreMatrix [][]int64, finalSeeds []*SeedDev, tempSeeds []*SeedDev, seedBuildHelper *seedHelper) []*SeedDev {
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
			seedBuildHelper.leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[seedBuildHelper.nodeIdx].SeqTwoBit, int(seedBuildHelper.nodePos), &read.Rainbow[seedBuildHelper.readOffset], readStart+seedBuildHelper.readOffset))
			tempSeeds = extendToTheRightDev(nodes[seedBuildHelper.nodeIdx], read, readStart-(seedBuildHelper.leftMatches-1), int(seedBuildHelper.nodePos)-(seedBuildHelper.leftMatches-1), true, tempSeeds)
			for _, seedBuildHelper.tempSeed = range tempSeeds {
				finalSeeds = append(finalSeeds, extendToTheLeftDev(nodes[seedBuildHelper.nodeIdx], read, seedBuildHelper.tempSeed)...)
			}
		}
		// do rev strand
		seedBuildHelper.seqKey = read.RainbowRc[seedBuildHelper.keyOffset].Seq[seedBuildHelper.keyIdx] >> seedBuildHelper.keyShift
		seedBuildHelper.currHits = seedHash[seedBuildHelper.seqKey]
		for _, seedBuildHelper.codedNodeCoord = range seedBuildHelper.currHits {
			seedBuildHelper.nodeIdx, seedBuildHelper.nodePos = numberToChromAndPos(seedBuildHelper.codedNodeCoord)
			seedBuildHelper.nodeOffset = int(seedBuildHelper.nodePos % basesPerInt)
			seedBuildHelper.readOffset = 31 - ((readStart - seedBuildHelper.nodeOffset + 31) % 32)

			seedBuildHelper.leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[seedBuildHelper.nodeIdx].SeqTwoBit, int(seedBuildHelper.nodePos), &read.RainbowRc[seedBuildHelper.readOffset], readStart+seedBuildHelper.readOffset))
			tempSeeds = extendToTheRightDev(nodes[seedBuildHelper.nodeIdx], read, readStart-(seedBuildHelper.leftMatches-1), int(seedBuildHelper.nodePos)-(seedBuildHelper.leftMatches-1), false, tempSeeds)
			finalSeeds = append(finalSeeds, tempSeeds...)
		}
	}
	if len(finalSeeds) > 100 {
		SortSeedDevByTotalLen(finalSeeds)
	} else {
		heapSortSeeds(finalSeeds)
	}
	return finalSeeds
}

func getLastPart(a *SeedDev) *SeedDev {
	for ; a.NextPart != nil; a = a.NextPart {
	}
	return a
}

func toTail(a *SeedDev) *SeedDev {
	for ; a.Next != nil; a = a.Next {
	}
	return a
}

func AlternateOrder(seeds []*SeedDev) []*SeedDev {
	var answer []*SeedDev = make([]*SeedDev, 0, len(seeds))
	for i, j := 0, len(seeds)-1; i < len(seeds)/2 || j > len(seeds)/2; i, j = i+1, j-1 {
		answer = append(answer, seeds[i])
		answer = append(answer, seeds[j])
	}
	return answer
}

func seedBedToSeed(a *SeedBed, currQPos uint32, posStrand bool) *SeedDev {
	if a == nil {
		return nil
	} else {
		return &SeedDev{TargetId: a.Id, TargetStart: a.Start, QueryStart: currQPos, Length: a.End - a.Start, PosStrand: posStrand, Next: seedBedToSeedDev(a.Next, currQPos+a.End-a.Start, posStrand)}
	}
}
