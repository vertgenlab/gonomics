package simpleGraph

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

// TODO: get rid of this when seedBed is eliminated
func seedBedToSeedDev(a *SeedBed, currQPos uint32, posStrand bool) *SeedDev {
	if a == nil {
		return nil
	} else {
		return &SeedDev{TargetId: a.Id, TargetStart: a.Start, QueryStart: currQPos, Length: a.End - a.Start, PosStrand: posStrand, Next: seedBedToSeedDev(a.Next, currQPos+a.End-a.Start, posStrand)}
	}
}

func extendToTheRight(node *Node, read *fastq.FastqBig, readStart int, nodeStart int, posStrand bool) []*SeedDev {
	const basesPerInt int = 32
	var nodeOffset int = nodeStart % basesPerInt
	var readOffset int = 31 - ((readStart - nodeOffset + 31) % 32)
	var rightMatches int = 0
	var currNode *SeedDev = nil
	var answer, nextParts []*SeedDev
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
			nextParts = extendToTheRight(node.Next[i].Dest, read, readStart+rightMatches, 0, posStrand)
			// if we aligned into the next node, make a seed for this node and point it to the next one
			for j = 0; j < len(nextParts); j++ {
				currNode = &SeedDev{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(rightMatches), PosStrand: posStrand, TotalLength: uint32(rightMatches) + nextParts[j].TotalLength, NextPart: nextParts[j], Next: nil}
				answer = append(answer, currNode)
			}
		}
	}

	// if the alignment did not go to another node, return the match for this node
	if len(answer) == 0 {
		currNode = &SeedDev{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(rightMatches), PosStrand: posStrand, TotalLength: uint32(rightMatches), NextPart: nil, Next: nil}
		answer = []*SeedDev{currNode}
	}

	return answer
}

func extendToTheLeft(node *Node, read *fastq.FastqBig, currPart *SeedDev) []*SeedDev {
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

func extendToTheLeftHelper(node *Node, read *fastq.FastqBig, nextPart *SeedDev) []*SeedDev {
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

func findSeedsInSmallMapWithMemPool(seedHash map[uint64][]uint64, nodes []*Node, read *fastq.FastqBig, seedLen int, perfectScore int64, scoreMatrix [][]int64) []*SeedDev {
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
			tempSeeds = extendToTheRight(nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), true)
			//log.Printf("After extendToTheRight fwd:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				finalSeeds = append(finalSeeds, extendToTheLeft(nodes[nodeIdx], read, tempSeed)...)
			}
			//log.Printf("After extendToTheLeft fwd:\n")
			//printSeedDev(finalSeeds)
			// TODO: Bring back speed optimizations once we are sure of correctness
			/*if leftMatches < 1 || rightMatches < seedLen {
				log.Fatalf("No matches found at seed location: %s %d, %d %d", dna.BasesToString(read.Seq), readStart, nodeIdx, nodePos)
			}*/

			/*if seedCouldBeBetter(int64(leftMatches+rightMatches-1), bestScore, perfectScore, int64(len(read.Seq)), 100, 90, -196, -296) {
				if poolHead != nil {
					currSeed = poolHead
					poolHead = poolHead.Next
				} else {
					currSeed = &SeedDev{}
				}
				currSeed.TargetId = uint32(nodeIdx)
				currSeed.TargetStart = uint32(int(nodePos) - leftMatches + 1)
				currSeed.QueryStart = uint32(readStart - leftMatches + 1)
				currSeed.Length = uint32(leftMatches + rightMatches - 1)
				currSeed.PosStrand = true
				currSeed.TotalLength = uint32(leftMatches + rightMatches - 1)
				currSeed.Next = hits
				hits = currSeed
				seedScore = scoreSeedFastqBig(currSeed, read, scoreMatrix)
				if seedScore > bestScore {
					bestScore = seedScore
				}
			}*/
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
			tempSeeds = extendToTheRight(nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), false)
			//log.Printf("After extendToTheRight rev:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				//log.Printf("tempSeed.QueryStart = %d\n", tempSeed.QueryStart)
				finalSeeds = append(finalSeeds, extendToTheLeft(nodes[nodeIdx], read, tempSeed)...)
			}
			//log.Printf("After extendToTheLeft rev:\n")
			//printSeedDev(finalSeeds)
			// TODO: bring back speed optimizations
			/*if leftMatches < 1 || rightMatches < seedLen {
				log.Fatalf("No matches found at seed location: %s %d, %d %d", dna.BasesToString(read.SeqRc), readStart, nodeIdx, nodePos)
			}*/
			/*if seedCouldBeBetter(int64(leftMatches+rightMatches-1), bestScore, perfectScore, int64(len(read.SeqRc)), 100, 90, -196, -296) {
				if poolHead != nil {
					currSeed = poolHead
					poolHead = poolHead.Next
				} else {
					currSeed = &SeedDev{}
				}
				currSeed.TargetId = uint32(nodeIdx)
				currSeed.TargetStart = uint32(int(nodePos) - leftMatches + 1)
				currSeed.QueryStart = uint32(readStart - leftMatches + 1)
				currSeed.Length = uint32(leftMatches + rightMatches - 1)
				currSeed.PosStrand = false
				currSeed.TotalLength = uint32(leftMatches + rightMatches - 1)
				currSeed.Next = hits
				hits = currSeed
				seedScore = scoreSeedFastqBig(currSeed, read, scoreMatrix)
				if seedScore > bestScore {
					bestScore = seedScore
				}
			}*/
		}
	}

	/*var finalHits, nextSeed *SeedDev
	for currSeed = hits; currSeed != nil; currSeed = nextSeed {
		nextSeed = currSeed.Next
		if seedCouldBeBetter(int64(currSeed.TotalLength), bestScore, perfectScore, int64(len(read.SeqRc)), 100, 90, -196, -296) {
			currSeed.Next = finalHits
			finalHits = currSeed
		} else {
			currSeed.Next = poolHead
			poolHead = currSeed
		}
	}

	*memoryPool = poolHead
	return hits*/
	//printSeedDev(finalSeeds)
	return finalSeeds
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

//TODO: get rid of this
func findSeedsInMapDev(seedHash map[uint64][]*SeedBed, read *fastq.Fastq, seedLen int, stepSize int, posStrand bool) []*SeedDev {
	var codedSeq uint64 = 0
	var hits []*SeedDev = make([]*SeedDev, 0)
	for subSeqStart := 0; subSeqStart < len(read.Seq)-seedLen+1; subSeqStart++ {
		if dna.CountBaseInterval(read.Seq, dna.N, subSeqStart, subSeqStart+seedLen) == 0 {
			codedSeq = dnaToNumber(read.Seq, subSeqStart, subSeqStart+seedLen)
			currHits := seedHash[codedSeq]
			for _, value := range currHits {
				hits = append(hits, seedBedToSeedDev(value, uint32(subSeqStart), posStrand))
			}
		}
	}
	//log.Printf("Total of %d hits.\n", len(hits))
	return hits
}

// need to handle neg strand
/*func findSeedsInMap(seedHash map[uint64][]*SeedBed, read *fastq.Fastq, seedLen int, stepSize int, posStrand bool) []*SeedDev {
	var codedSeq uint64 = 0
	var prevHits []*SeedDev = make([]*SeedDev, 0)
	var allHits []*SeedDev = make([]*SeedDev, 0)
	for initOffset := 0; initOffset < stepSize; initOffset++ {
		for subSeqStart := initOffset; subSeqStart < len(read.Seq)-seedLen+1; subSeqStart += stepSize {
			if dna.CountBaseInterval(read.Seq, dna.N, subSeqStart, subSeqStart+seedLen) == 0 {
				codedSeq = dnaToNumber(read.Seq, subSeqStart, subSeqStart+seedLen)
				//fmt.Printf("Coded seq is:%d, seedLength:%d\n", codedSeq, seedLen)
				currHits := seedHash[codedSeq]
				//log.Printf("At position %d, we found %d hits.\n", subSeqStart, len(currHits))
				noMerge, merged := mergeSeedLists(prevHits, currHits, uint32(subSeqStart), posStrand)
				allHits = append(allHits, noMerge...)
				prevHits = merged
			}
		}
		allHits = append(allHits, prevHits...)
		prevHits = prevHits[0:0]
	}
	//log.Printf("Total of %d hits.\n", len(allHits))
	return allHits
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
