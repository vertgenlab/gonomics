package genomeGraph

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/numbers"
)

// extendToTheRight extends a seed to the right, considering possible continuations in connected nodes.
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
	if currPart.QueryStart == 0 || currPart.TargetStart > 0 {
		return []*SeedDev{currPart}
	}

	answer := make([]*SeedDev, 0)
	for _, prev := range node.Prev {
		var readBase dna.Base
		if currPart.PosStrand {
			readBase = dnaTwoBit.GetBase(&read.Rainbow[0], uint(currPart.QueryStart)-1)
		} else {
			readBase = dnaTwoBit.GetBase(&read.RainbowRc[0], uint(currPart.QueryStart)-1)
		}
		if readBase == dnaTwoBit.GetBase(prev.Dest.SeqTwoBit, uint(prev.Dest.SeqTwoBit.Len)-1) {
			prevParts := extendToTheLeftHelper(prev.Dest, read, currPart)
			answer = append(answer, prevParts...)
		}
	}

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
		answer := make([]*SeedDev, 0)
		for _, prev := range node.Prev {
			prevParts := extendToTheLeftHelper(prev.Dest, read, currPart)
			answer = append(answer, prevParts...)
		}
		return answer
	}

	return []*SeedDev{currPart}
}

func findSeedsInSmallMapWithMemPool(seedHash map[uint64][]uint64, nodes []Node, read fastq.FastqBig, seedLen int, perfectScore int64, scoreMatrix [][]int64) []*SeedDev {
	// Pre-compute constants
	const basesPerInt = 32
	keyShift := uint(64 - (seedLen * 2))

	// Initialize slice to store final seeds
	finalSeeds := make([]*SeedDev, 0)

	// Loop through each start position in the read
	for readStart := 0; readStart <= len(read.Seq)-seedLen; readStart++ {
		// Forward strand processing
		processStrand(seedHash, nodes, read, seedLen, &finalSeeds, readStart, true, keyShift)

		// Reverse strand processing
		processStrand(seedHash, nodes, read, seedLen, &finalSeeds, readStart, false, keyShift)
	}

	return finalSeeds
}

func processStrand(seedHash map[uint64][]uint64, nodes []Node, read fastq.FastqBig, seedLen int, finalSeeds *[]*SeedDev, readStart int, posStrand bool, keyShift uint) {
	keyIdx := (readStart + 31) / 32
	keyOffset := 31 - ((readStart + 31) % 32)
	var seqKey uint64
	if posStrand {
		seqKey = read.Rainbow[keyOffset].Seq[keyIdx] >> keyShift
	} else {
		seqKey = read.RainbowRc[keyOffset].Seq[keyIdx] >> keyShift
	}

	currHits := seedHash[seqKey]
	for _, codedNodeCoord := range currHits {
		nodeIdx, nodePos := numberToChromAndPos(codedNodeCoord)
		tempSeeds := extendToTheRight(&nodes[nodeIdx], read, readStart, int(nodePos), posStrand)
		for _, tempSeed := range tempSeeds {
			*finalSeeds = append(*finalSeeds, extendToTheLeft(&nodes[nodeIdx], read, tempSeed)...)
		}
	}
}

// does not support extending along edges, so commented out
/*func findSeedsInSmallMap(seedHash map[uint64][]uint64, nodes []*Node, read *fastq.FastqBig, seedLen int, perfectScore int64, scoreMatrix [][]int64) []*SeedDev {
	var currHits []uint64
	var leftMatches int = 0
	var nodeIdx int64 = 0
	var readStart, keyIdx, keyOffset int
	var seqKey, codedNodeCoord uint64
	var keyShift uint
	var nodePos int64
	var nodeOffset, basesPerInt, readOffset int
	var tempSeeds []*SeedDev
	var tempSeed *SeedDev
	var finalSeeds []*SeedDev

	log.Printf("About to find seeds:\n")
	for readStart = 0; readStart < len(read.Seq)-seedLen+1; readStart++ {
		keyIdx = (readStart + 31) / 32
		keyOffset = 31 - ((readStart + 31) % 32)

		// do fwd strand
		seqKey = read.Rainbow[keyOffset].Seq[keyIdx] >> keyShift
		currHits = seedHash[seqKey]
		log.Printf("Found %d forward hits\n", len(currHits))
		for _, codedNodeCoord = range currHits {
			nodeIdx, nodePos = numberToChromAndPos(codedNodeCoord)
			nodeOffset = int(nodePos) % basesPerInt
			readOffset = 31 - ((readStart - nodeOffset + 31) % 32)

			leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.Rainbow[readOffset], readStart+readOffset))
			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.Rainbow[readOffset], readStart+readOffset)
			tempSeeds = extendToTheRight(nodes[nodeIdx], read, readStart-leftMatches, int(nodePos)-leftMatches, true)
			log.Printf("After extendToTheRight fwd:\n")
			printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				finalSeeds = append(finalSeeds, extendToTheLeft(nodes[nodeIdx], read, tempSeed)...)
			}
			log.Printf("After extendToTheLeft fwd:\n")
			printSeedDev(finalSeeds)
		}

		// do rev strand
		seqKey = read.RainbowRc[keyOffset].Seq[keyIdx] >> keyShift
		currHits = seedHash[seqKey]
		log.Printf("Found %d reverse hits\n", len(currHits))
		for _, codedNodeCoord = range currHits {
			nodeIdx, nodePos = numberToChromAndPos(codedNodeCoord)
			nodeOffset = int(nodePos) % basesPerInt
			readOffset = 31 - ((readStart - nodeOffset + 31) % 32)

			leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.RainbowRc[readOffset], readStart+readOffset))
			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.RainbowRc[readOffset], readStart+readOffset)
			tempSeeds = extendToTheRight(nodes[nodeIdx], read, readStart-leftMatches, int(nodePos)-leftMatches, false)
			log.Printf("After extendToTheRight rev:\n")
			printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				finalSeeds = append(finalSeeds, extendToTheLeft(nodes[nodeIdx], read, tempSeed)...)
			}
			log.Printf("After extendToTheLeft Rev:\n")
			printSeedDev(finalSeeds)
		}
	}
	return finalSeeds
	// TODO: Bring back speed optimizations
	var badIdx = len(hits) - 1
	var badCount = 0
	for i := 0; i <= badIdx; {
		if !seedCouldBeBetter(int64(hits[i].TotalLength), bestScore, perfectScore, int64(len(read.SeqRc)), 100, 90, -196, -296) {
			hits[i], hits[badIdx] = hits[badIdx], hits[i]
			badIdx--
			badCount++
		} else {
			i++
		}
	}
	return hits[0:(badIdx + 1)]
}*/

/*func findSeedsInSlice(seedHash [][]uint64, read *fastq.Fastq, seedLen int, posStrand bool) []*SeedDev {
	var codedSeq uint64 = 0
	var allHits []*SeedDev = make([]*SeedDev, 0)
	for subSeqStart := 0; subSeqStart < len(read.Seq)-seedLen+1; subSeqStart++ {
		if dna.CountBaseInterval(read.Seq, dna.N, subSeqStart, subSeqStart+seedLen) == 0 {
			codedSeq = dnaToNumber(read.Seq, subSeqStart, subSeqStart+seedLen)
			//fmt.Printf("Coded seq is:%d, seedLength:%d\n", codedSeq, seedLen)
			currHits := seedHash[codedSeq]
			noMerge, merged := mergeSeedLists(prevHits, currHits, uint32(subSeqStart), posStrand)
			allHits = append(allHits, noMerge...)
			prevHits = merged
		}

	}
	allHits = append(allHits, prevHits...)

	return allHits
}*/
