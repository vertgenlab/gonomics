package genomeGraph

import (
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
		for _, next := range node.Next {
			nextParts := extendToTheRight(next.Dest, read, readStart+rightMatches, 0, posStrand)
			for _, part := range nextParts {
				currNode := &SeedDev{
					TargetId:    node.Id,
					TargetStart: uint32(nodeStart),
					QueryStart:  uint32(readStart),
					Length:      uint32(rightMatches),
					PosStrand:   posStrand,
					TotalLength: uint32(rightMatches) + part.TotalLength,
					NextPart:    part,
				}
				answer = append(answer, currNode)
			}
		}
	}

	if len(answer) == 0 {
		answer = append(answer, &SeedDev{
			TargetId:    node.Id,
			TargetStart: uint32(nodeStart),
			QueryStart:  uint32(readStart),
			Length:      uint32(rightMatches),
			PosStrand:   posStrand,
			TotalLength: uint32(rightMatches),
		})
	}

	return answer
}

// extendToTheLeft extends a seed to the left, considering possible continuations in connected nodes.
func extendToTheLeft(node *Node, read fastq.FastqBig, currPart *SeedDev) []*SeedDev {
	// Immediate return if extension to the left is not possible
	if currPart.QueryStart == 0 || currPart.TargetStart > 0 {
		return []*SeedDev{currPart}
	}

	// Precompute read bases for both strands
	readBaseFwd := dnaTwoBit.GetBase(&read.Rainbow[0], uint(currPart.QueryStart)-1)
	readBaseRev := dnaTwoBit.GetBase(&read.RainbowRc[0], uint(currPart.QueryStart)-1)

	var answer []*SeedDev
	for _, prev := range node.Prev {
		prevNodeLastBase := dnaTwoBit.GetBase(prev.Dest.SeqTwoBit, uint(prev.Dest.SeqTwoBit.Len)-1)

		// Check if the last base of the previous node matches the base before the current part's start
		if (currPart.PosStrand && readBaseFwd == prevNodeLastBase) || (!currPart.PosStrand && readBaseRev == prevNodeLastBase) {
			// Extend seed to the left within the previous node
			prevParts := extendToTheLeftHelper(prev.Dest, read, currPart)
			answer = append(answer, prevParts...)
		}
	}

	// Return the current part if no extension is possible
	if len(answer) == 0 {
		return []*SeedDev{currPart}
	}
	return answer
}

// extendToTheLeftHelper assists in extending a seed to the left within a specific node.
func extendToTheLeftHelper(node *Node, read fastq.FastqBig, nextPart *SeedDev) []*SeedDev {
	nodePos := node.SeqTwoBit.Len - 1
	readPos := int(nextPart.QueryStart) - 1
	nodeOffset := nodePos % 32
	readOffset := 31 - ((readPos - nodeOffset + 31) % 32)
	var leftMatches int
	if nextPart.PosStrand {
		leftMatches = numbers.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.Rainbow[readOffset], readPos+readOffset))
	} else {
		leftMatches = numbers.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.RainbowRc[readOffset], readPos+readOffset))
	}

	currPart := &SeedDev{
		TargetId:    node.Id,
		TargetStart: uint32(nodePos - (leftMatches - 1)),
		QueryStart:  uint32(readPos - (leftMatches - 1)),
		Length:      uint32(leftMatches),
		PosStrand:   nextPart.PosStrand,
		TotalLength: uint32(leftMatches) + nextPart.TotalLength,
		NextPart:    nextPart,
	}

	if currPart.QueryStart > 0 && currPart.TargetStart == 0 {
		answer := make([]*SeedDev, 0, len(node.Prev))
		for _, prev := range node.Prev {
			answer = append(answer, extendToTheLeftHelper(prev.Dest, read, currPart)...)
		}
		return answer
	}

	return []*SeedDev{currPart}
}

func findSeedsInSmallMapWithMemPool(seedHash map[uint64][]uint64, nodes []Node, read fastq.FastqBig, seedLen int, perfectScore int64, scoreMatrix [][]int64) []*SeedDev {
	keyShift := uint(64 - (seedLen * 2))

	// Initialize slice to store final seeds
	finalSeeds := make([]*SeedDev, 0)
	// Loop through each start position in the read
	for readStart := 0; readStart <= len(read.Seq)-seedLen; readStart++ {
		// Forward strand processing
		processStrand(seedHash, nodes, read, seedLen, &finalSeeds, readStart, keyShift)
		// Reverse strand processing
	}

	return finalSeeds
}

func processStrand(seedHash map[uint64][]uint64, nodes []Node, read fastq.FastqBig, seedLen int, finalSeeds *[]*SeedDev, readStart int, keyShift uint) {
	keyIdx := (readStart + 31) / 32
	keyOffset := 31 - ((readStart + 31) % 32)

	var fwdSeqKey, revSeqKey uint64 = read.Rainbow[keyOffset].Seq[keyIdx] >> keyShift, read.RainbowRc[keyOffset].Seq[keyIdx] >> keyShift
	var nodeIdx, nodePos int64

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
