package simpleGraph

import (
	"github.com/vertgenlab/gonomics/common"
	//"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	//"log"
)

/*func seedBedToSeedDev(a *SeedBed, currQPos uint32, posStrand bool) *SeedDev {
	if a == nil {
		return nil
	} else {
		return &SeedDev{TargetId: a.Id, TargetStart: a.Start, QueryStart: currQPos, Length: a.End - a.Start, PosStrand: posStrand, Next: seedBedToSeedDev(a.Next, currQPos+a.End-a.Start, posStrand)}
	}
}*/

/*func mergeIntoSeed(aTail *SeedDev, b *SeedBed, currQPos uint32, posStrand bool) bool {
	if aTail.TargetId == b.Id &&
		common.MaxUint32(aTail.TargetStart, b.Start) <= common.MinUint32(aTail.TargetStart+aTail.Length, b.End) &&
		common.MaxUint32(aTail.QueryStart, currQPos) <= common.MinUint32(aTail.QueryStart+aTail.Length, currQPos+b.End-b.Start) &&
		b.Start > aTail.TargetStart && currQPos > aTail.QueryStart && (b.Start-aTail.TargetStart == currQPos-aTail.QueryStart) {
		extension := b.End - (aTail.TargetStart + aTail.Length)
		aTail.Length += extension
		if aTail.Next != nil {
			//I think this may happen when a seed can be extended along two different edges.  The preceeding seed nodes would need to be duplicated.
			log.Fatal("Error: Trouble merging alignment seeds, Next should be nil.\n")
		}
		aTail.Next = seedBedToSeedDev(b.Next, currQPos+b.End-b.Start, posStrand)
		return true
	} else {
		if b.Next == nil {
			return false
		} else {
			return mergeIntoSeed(aTail, b.Next, currQPos+b.End-b.Start, posStrand)
		}
	}
}*/

/*func getTails(a []*SeedDev) []*SeedDev {
	tails := make([]*SeedDev, len(a))
	var tail *SeedDev = nil
	for i := 0; i < len(a); i++ {
		for tail = a[0]; tail.Next != nil; tail = tail.Next {
		}
		tails[i] = tail
	}
	return tails
}*/

/*func mergeSeedLists(lastPosition []*SeedDev, currPosition []*SeedBed, currQPos uint32, posStrand bool) ([]*SeedDev, []*SeedDev) {
	noMerge := make([]*SeedDev, 0)
	merged := make([]*SeedDev, 0)
	var lastPositionTail *SeedDev = nil
	var usedLastPosition bool = false
	usedCurrent := make([]bool, len(currPosition))
	//log.Printf("len(lastPosition)=%d, len(currPosition)=%d\n", len(lastPosition), len(currPosition))
	for j := 0; j < len(lastPosition); j++ {
		for lastPositionTail = lastPosition[j]; lastPositionTail.Next != nil; lastPositionTail = lastPositionTail.Next {
		}
		usedLastPosition = false
		for k := 0; k < len(currPosition); k++ {
			if mergeIntoSeed(lastPositionTail, currPosition[k], currQPos, posStrand) {
				usedLastPosition = true
				usedCurrent[k] = true
			}
		}
		if usedLastPosition {
			merged = append(merged, lastPosition[j])
		} else {
			noMerge = append(noMerge, lastPosition[j])
		}
	}
	for k := 0; k < len(usedCurrent); k++ {
		if !usedCurrent[k] {
			merged = append(merged, seedBedToSeedDev(currPosition[k], currQPos, posStrand))
		}
	}
	//log.Printf("len(noMerge)=%d, len(merged)=%d\n", len(noMerge), len(merged))
	//log.Printf("noMerge\n")
	//printSeedDev(noMerge)
	//log.Printf("merged\n")
	//printSeedDev(merged)
	return noMerge, merged
}*/

func findSeedsInSmallMapWithMemPool(seedHash map[uint64][]uint64, nodes []*Node, read *fastq.FastqBig, seedLen int, perfectScore int64, memoryPool **SeedDev) *SeedDev {
	var hits *SeedDev
	var currHits []uint64
	var codedPos uint64
	var currSeed *SeedDev
	var leftMatches, rightMatches int = 0, 0
	var idx, offset int = 0, 0
	var nodeIdx, pos int64 = 0, 0
	var bestScore, seedScore int64 = 0, 0
	var poolHead *SeedDev = *memoryPool

	for readStart := 0; readStart < len(read.Seq)-seedLen+1; readStart++ {
		idx = (readStart + 31) / 32
		offset = 31 - ((readStart + 31) % 32)

		// do fwd strand
		currHits = seedHash[read.Rainbow[offset].Seq[idx]]
		for _, codedPos = range currHits {
			nodeIdx, pos = numberToChromAndPos(codedPos)
			leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(pos), read.Rainbow[offset], readStart+offset))
			rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(pos), read.Rainbow[offset], readStart+offset)
			if seedCouldBeBetter(int64(leftMatches+rightMatches-1), bestScore, perfectScore, int64(len(read.Seq)), 100, 90, -196, -296) {
				if poolHead != nil {
					currSeed = poolHead
					poolHead = poolHead.Next
				} else {
					currSeed = &SeedDev{}
				}
				currSeed.TargetId = uint32(nodeIdx)
				currSeed.TargetStart = uint32(int(pos) - leftMatches + 1)
				currSeed.QueryStart = uint32(readStart - leftMatches + 1)
				currSeed.Length = uint32(leftMatches + rightMatches - 1)
				currSeed.PosStrand = true
				currSeed.TotalLength = uint32(leftMatches + rightMatches - 1)
				currSeed.Next = hits
				hits = currSeed
				seedScore = scoreSeedFastqBig(currSeed, read)
				if seedScore > bestScore {
					bestScore = seedScore
				}
			}
		}

		// do rev strand
		currHits = seedHash[read.RainbowRc[offset].Seq[idx]]
		for _, codedPos = range currHits {
			nodeIdx, pos = numberToChromAndPos(codedPos)
			leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(pos), read.RainbowRc[offset], readStart+offset))
			rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(pos), read.RainbowRc[offset], readStart+offset)
			if seedCouldBeBetter(int64(leftMatches+rightMatches-1), bestScore, perfectScore, int64(len(read.SeqRc)), 100, 90, -196, -296) {
				if poolHead != nil {
					currSeed = poolHead
					poolHead = poolHead.Next
				} else {
					currSeed = &SeedDev{}
				}
				currSeed.TargetId = uint32(nodeIdx)
				currSeed.TargetStart = uint32(int(pos) - leftMatches + 1)
				currSeed.QueryStart = uint32(readStart - leftMatches + 1)
				currSeed.Length = uint32(leftMatches + rightMatches - 1)
				currSeed.PosStrand = false
				currSeed.TotalLength = uint32(leftMatches + rightMatches - 1)
				currSeed.Next = hits
				hits = currSeed
				seedScore = scoreSeedFastqBig(currSeed, read)
				if seedScore > bestScore {
					bestScore = seedScore
				}
			}
		}
	}

	var finalHits, nextSeed *SeedDev
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
	return hits
}

func findSeedsInSmallMap(seedHash map[uint64][]uint64, nodes []*Node, read *fastq.FastqBig, seedLen int, perfectScore int64) []*SeedDev {
	var hits []*SeedDev = make([]*SeedDev, 0)
	var currHits []uint64
	var codedPos uint64
	var currSeed *SeedDev
	var leftMatches, rightMatches int = 0, 0
	var idx, offset int = 0, 0
	var nodeIdx, pos int64 = 0, 0
	var bestScore, seedScore int64 = 0, 0

	for readStart := 0; readStart < len(read.Seq)-seedLen+1; readStart++ {
		idx = (readStart + 31) / 32
		offset = 31 - ((readStart + 31) % 32)

		// do fwd strand
		currHits = seedHash[read.Rainbow[offset].Seq[idx]]
		for _, codedPos = range currHits {
			nodeIdx, pos = numberToChromAndPos(codedPos)
			leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(pos), read.Rainbow[offset], readStart+offset))
			rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(pos), read.Rainbow[offset], readStart+offset)
			if seedCouldBeBetter(int64(leftMatches+rightMatches-1), bestScore, perfectScore, int64(len(read.Seq)), 100, 90, -196, -296) {
				currSeed = &SeedDev{TargetId: uint32(nodeIdx), TargetStart: uint32(int(pos) - leftMatches + 1), QueryStart: uint32(readStart - leftMatches + 1), Length: uint32(leftMatches + rightMatches - 1), PosStrand: true, TotalLength: uint32(leftMatches + rightMatches - 1), Next: nil}
				seedScore = scoreSeedFastqBig(currSeed, read)
				if seedScore > bestScore {
					bestScore = seedScore
				}
				hits = append(hits, currSeed)
			}
		}

		// do rev strand
		currHits = seedHash[read.RainbowRc[offset].Seq[idx]]
		for _, codedPos = range currHits {
			nodeIdx, pos = numberToChromAndPos(codedPos)
			leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(pos), read.RainbowRc[offset], readStart+offset))
			rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(pos), read.RainbowRc[offset], readStart+offset)
			if seedCouldBeBetter(int64(leftMatches+rightMatches-1), bestScore, perfectScore, int64(len(read.SeqRc)), 100, 90, -196, -296) {
				currSeed = &SeedDev{TargetId: uint32(nodeIdx), TargetStart: uint32(int(pos) - leftMatches + 1), QueryStart: uint32(readStart - leftMatches + 1), Length: uint32(leftMatches + rightMatches - 1), PosStrand: false, TotalLength: uint32(leftMatches + rightMatches - 1), Next: nil}
				seedScore = scoreSeedFastqBig(currSeed, read)
				if seedScore > bestScore {
					bestScore = seedScore
				}
				hits = append(hits, currSeed)
			}
		}
	}
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
}

/*func findSeedsInMapDev(seedHash map[uint64][]*SeedBed, read *fastq.Fastq, seedLen int, stepSize int, posStrand bool) []*SeedDev {
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
}*/

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

// need to handle neg strand
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
