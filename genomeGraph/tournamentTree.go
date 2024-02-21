package genomeGraph

import (
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/numbers"
)

const basesPerInt = 32

func extendToTheRight(node *Node, read fastq.FastqBig, readStart, nodeStart int, posStrand bool) []*SeedDev {
	const basesPerInt = 32
	nodeOffset := nodeStart % basesPerInt
	readOffset := 31 - ((readStart - nodeOffset + 31) % 32)
	var rightMatches int
	if posStrand {
		rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, &read.Rainbow[readOffset], readStart+readOffset)
	} else {
		rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, &read.RainbowRc[readOffset], readStart+readOffset)
	}

	if rightMatches == 0 {
		return nil
	}

	answer := make([]*SeedDev, 0)
	if readStart+rightMatches < len(read.Seq) && nodeStart+rightMatches == node.SeqTwoBit.Len && len(node.Next) > 0 {
		answer = appendSeedExtensions(node, read, readStart, nodeStart, posStrand, rightMatches, answer)
	}
	if len(answer) == 0 {
		answer = append(answer, createSeed(node.Id, readStart, nodeStart, rightMatches, posStrand))
	}
	return answer
}

func generateLeftExtensionSeeds(node *Node, read fastq.FastqBig, currPart *SeedDev) []*SeedDev {
	seeds := make([]*SeedDev, 0, len(node.Prev))
	for _, prev := range node.Prev {
		seeds = append(seeds, extendToTheLeftHelper(prev.Dest, read, currPart)...)
	}
	return seeds
}

func extendToTheLeft(node *Node, read fastq.FastqBig, currPart *SeedDev) []*SeedDev {
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
				prevParts = extendToTheLeftHelper(node.Prev[i].Dest, read, currPart)
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

func extendToTheLeftHelper(node *Node, read fastq.FastqBig, nextPart *SeedDev) []*SeedDev {
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
		leftMatches = numbers.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.Rainbow[readOffset], readPos+readOffset))
	} else {
		leftMatches = numbers.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.RainbowRc[readOffset], readPos+readOffset))
	}

	if leftMatches == 0 {
		log.Fatal("Error: should not have zero matches to the left\n")
	}

	currPart = &SeedDev{TargetId: node.Id, TargetStart: uint32(nodePos - (leftMatches - 1)), QueryStart: uint32(readPos - (leftMatches - 1)), Length: uint16(leftMatches), PosStrand: nextPart.PosStrand, TotalLength: uint16(leftMatches) + nextPart.TotalLength, NextPart: nextPart}

	// we went all the way to end and there might be more
	if currPart.QueryStart > 0 && currPart.TargetStart == 0 {
		for i = 0; i < len(node.Prev); i++ {
			if nextPart.PosStrand {
				readBase = dnaTwoBit.GetBase(&read.Rainbow[0], uint(currPart.QueryStart)-1)
			} else {
				readBase = dnaTwoBit.GetBase(&read.RainbowRc[0], uint(currPart.QueryStart)-1)
			}
			if readBase == dnaTwoBit.GetBase(node.Prev[i].Dest.SeqTwoBit, uint(node.Prev[i].Dest.SeqTwoBit.Len)-1) {
				prevParts = extendToTheLeftHelper(node.Prev[i].Dest, read, currPart)
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

func findSeedsInSmallMapWithMemPool(seedHash map[uint64][]uint64, nodes []Node, read fastq.FastqBig, seedLen int, perfectScore int64, scoreMatrix [][]int64) []*SeedDev {
	const basesPerInt int64 = 32
	var currHits []uint64
	var codedNodeCoord uint64
	var leftMatches int = 0
	var keyIdx, keyOffset, readOffset, nodeOffset int = 0, 0, 0, 0
	var nodeIdx, nodePos int64 = 0, 0
	//var poolHead *SeedDev = *memoryPool
	var seqKey uint64
	var keyShift uint = 64 - (uint(seedLen) * 2)
	var tempSeeds, finalSeeds []*SeedDev
	var tempSeed *SeedDev

	for readStart := 0; readStart < len(read.Seq)-seedLen+1; readStart++ {
		keyIdx = (readStart + 31) / 32
		keyOffset = 31 - ((readStart + 31) % 32)

		// do fwd strand
		seqKey = read.Rainbow[keyOffset].Seq[keyIdx] >> keyShift
		currHits = seedHash[seqKey]
		/*if len(currHits) > 0 {
			log.Printf(" %d hits\n", len(currHits))
		}*/
		for _, codedNodeCoord = range currHits {
			nodeIdx, nodePos = numberToChromAndPos(codedNodeCoord)
			nodeOffset = int(nodePos % basesPerInt)
			readOffset = 31 - ((readStart - nodeOffset + 31) % 32)

			leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), &read.Rainbow[readOffset], readStart+readOffset))
			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.Rainbow[readOffset], readStart+readOffset)
			tempSeeds = extendToTheRight(&nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), true)
			//log.Printf("After extendToTheRight fwd:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				finalSeeds = append(finalSeeds, extendToTheLeft(&nodes[nodeIdx], read, tempSeed)...)
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

			leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), &read.RainbowRc[readOffset], readStart+readOffset))
			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.RainbowRc[readOffset], readStart+readOffset)
			tempSeeds = extendToTheRight(&nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), false)
			//log.Printf("After extendToTheRight rev:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				//log.Printf("tempSeed.QueryStart = %d\n", tempSeed.QueryStart)
				finalSeeds = append(finalSeeds, extendToTheLeft(&nodes[nodeIdx], read, tempSeed)...)
			}
		}
	}

	return finalSeeds
}

// generateSeeds generates the seeds for extend to the right operations.
func generateSeeds(node *Node, read fastq.FastqBig, readStart, nodeStart int, posStrand bool, rightMatches int) []*SeedDev {
	seeds := make([]*SeedDev, 0)
	isEndOfNode := nodeStart+rightMatches == node.SeqTwoBit.Len

	if isEndOfNode && len(node.Next) > 0 {
		seeds = appendSeedExtensions(node, read, readStart, nodeStart, posStrand, rightMatches, seeds)
	}

	if len(seeds) == 0 {
		seeds = append(seeds, createSeed(node.Id, readStart, nodeStart, rightMatches, posStrand))
	}

	return seeds
}

// appendSeedExtensions appends the seeds of the extensions constructed to the right of the current one.
func appendSeedExtensions(node *Node, read fastq.FastqBig, readStart, nodeStart int, posStrand bool, rightMatches int, seeds []*SeedDev) []*SeedDev {
	for _, next := range node.Next {
		nextParts := extendToTheRight(next.Dest, read, readStart+rightMatches, 0, posStrand)
		for _, part := range nextParts {
			currNode := createSeed(node.Id, readStart, nodeStart, rightMatches, posStrand)
			currNode.TotalLength += part.TotalLength
			currNode.NextPart = part
			seeds = append(seeds, currNode)
		}
	}
	return seeds
}

func processStrand(seedHash map[uint64][]uint64, nodes []Node, read fastq.FastqBig, seedLen int, finalSeeds *[]*SeedDev, readStart int, keyShift uint) {
	keyIdx := (readStart + 31) / 32
	keyOffset := 31 - ((readStart + 31) % 32)

	var fwdSeqKey, revSeqKey uint64 = read.Rainbow[keyOffset].Seq[keyIdx] >> keyShift, read.RainbowRc[keyOffset].Seq[keyIdx] >> keyShift
	var nodeIdx, nodePos int64

	// Forward + Reverse strand processing
	fwdSeedHits := seedHash[fwdSeqKey]
	revSeedHits := seedHash[revSeqKey]

	fwdSize, revSize := len(fwdSeedHits), len(revSeedHits)

	for index := 0; index < fwdSize || index < revSize; index++ {
		if index < fwdSize {
			nodeIdx, nodePos = numberToChromAndPos(fwdSeedHits[index])
			processSeed(&nodes[nodeIdx], read, readStart, int(nodePos), true, finalSeeds)
		}
		if index < revSize {
			nodeIdx, nodePos = numberToChromAndPos(revSeedHits[index])
			processSeed(&nodes[nodeIdx], read, readStart, int(nodePos), false, finalSeeds)
		}
	}
}

func processSeed(node *Node, read fastq.FastqBig, readStart int, nodePos int, posStrand bool, finalSeeds *[]*SeedDev) {
	tempSeeds := extendToTheRight(node, read, readStart, nodePos, posStrand)
	for _, tempSeed := range tempSeeds {
		*finalSeeds = append(*finalSeeds, extendToTheLeft(node, read, tempSeed)...)
	}
}
