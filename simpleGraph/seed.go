package simpleGraph

import (
	//"github.com/vertgenlab/gonomics/common"
	//"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
	"sort"
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

func GraphDictionary(seeds []*SeedDev, gg *SimpleGraph, read *fastq.Fastq) []*SeedDev {
	var hash []*SeedDev
	for i := 0; i < len(seeds); i++ {
		hash = append(hash, extendSeedTogether(seeds[i], gg, read)...)
	}
	//SortBlastz(hash, read)
	//if len(hash) > 248 {
		//log.Printf("Length of hash for read:%s is %d\n", read.Name, len(hash))
	//	hash = hash[:30]
	//}
	return hash
}
/*
// TODO: this does not take into account breaking up seeds by gaps instead of mismatches
// similar calculations could also be used as the parameters to a banded alignment
func isSeedBetter(currIndex int, hash []*SeedDev, currBestScore int64, perfectScore int64, queryLen int64, maxMatch int64, minMatch int64, leastSevereMismatch int64, leastSevereMatchMismatchChange int64) bool {
	seedLen := int64(sumLen(hash[currIndex]))
	seeds := queryLen / (seedLen + 1)
	remainder := queryLen % (seedLen + 1)
	if currIndex > 248 {
		//log.Printf("seed along with whole seeds, but not the remainder: %d\n", seedLen*seeds*maxMatch+seeds*leastSevereMismatch)
		return false
	} else if currIndex > 30 && currBestScore > 11000 {
		//log.Printf("Index where we found a score of 10000: %d\n", currIndex)
		return false
	} else {
		// seed by itself could be best
		if seedLen*maxMatch >= currBestScore &&
			perfectScore-((queryLen-seedLen)*minMatch) >= currBestScore {
			return true
			// seed along with whole seeds, but not the remainder
		} else if seedLen*seeds*maxMatch+seeds*leastSevereMismatch >= currBestScore &&
			perfectScore-remainder*minMatch+seeds*leastSevereMatchMismatchChange >= currBestScore {
			return true
			// seed along with whole seeds, as well as both remainders
		} else if seedLen*seeds*maxMatch+remainder*maxMatch+(seeds+1)*leastSevereMismatch >= currBestScore &&
			perfectScore+(seeds+1)*leastSevereMatchMismatchChange >= currBestScore {
			return true
		} else {
			return false
		}
	}
}*/

func extendCurrSeed(seed *SeedDev, gg *SimpleGraph, read *fastq.Fastq, left bool, right bool) {
	var newTStart, newQStart, newTEnd, newQEnd int32 = int32(seed.TargetStart) - 1, int32(seed.QueryStart) - 1, int32(seed.TargetStart + seed.Length), int32(seed.QueryStart + seed.Length)
	//check to see if begining is at index zero, if so do something like SeedDev.Prev
	//if newStart < 0
	if left {
		for ; newTStart >= 0 && newQStart >= 0 && (gg.Nodes[seed.TargetId].Seq[newTStart] == read.Seq[newQStart]); newTStart, newQStart = newTStart-1, newQStart-1 {
			seed.TargetStart = uint32(newTStart)
			seed.QueryStart = uint32(newQStart)
			seed.Length++
		}
	}
	if right {
		for ; int(newTEnd) < len(gg.Nodes[seed.TargetId].Seq) && int(newQEnd) < len(read.Seq) && (gg.Nodes[seed.TargetId].Seq[newTEnd] == read.Seq[newQEnd]); newTEnd, newQEnd = newTEnd+1, newQEnd+1 {
			seed.Length++
		}
	}
}

func toTheRight(seed *SeedDev, gg *SimpleGraph, read *fastq.Fastq) []*SeedDev {
	//log.Printf("Depth of call is: %d for seed: %d", depth, seed.TargetId)
	var answer []*SeedDev
	extendCurrSeed(seed, gg, read, false, true)
	var newTEnd, newQEnd int32 = int32(seed.TargetStart + seed.Length), int32(seed.QueryStart + seed.Length)
	if (int(newTEnd) >= len(gg.Nodes[seed.TargetId].Seq) && int(newQEnd) < len(read.Seq)) && (len(gg.Nodes[seed.TargetId].Next) > 0) {
		var newTStart int32 = 0
		var newQStart int32 = newQEnd
		var edgeSeeds []*SeedDev
		var e int
		for _, next := range gg.Nodes[seed.TargetId].Next {
			//log.Printf("Number of nodes to check %d\n", len(gg.Nodes[seed.TargetId].Next))
			if len(next.Dest.Seq) > 0 {
				if next.Dest.Seq[newTStart] == read.Seq[newQStart] {
					nextSeed := &SeedDev{TargetId: next.Dest.Id, TargetStart: uint32(newTStart), QueryStart: uint32(newQStart), Length: 1, PosStrand: seed.PosStrand, Next: nil, Prev: nil}
					edgeSeeds = toTheRight(nextSeed, gg, read)
					for e = 0; e < len(edgeSeeds); e++ {
						currSeed := &SeedDev{TargetId: seed.TargetId, TargetStart: seed.TargetStart, QueryStart: seed.QueryStart, Length: seed.Length, PosStrand: seed.PosStrand, Next: nil, Prev: nil}
						currSeed.Next = edgeSeeds[e]
						answer = append(answer, currSeed)
					}
				}
			}
		}
	} else {
		answer = append(answer, seed)
	}
	return answer
}

func toTheLeft(seed *SeedDev, gg *SimpleGraph, read *fastq.Fastq) []*SeedDev {
	var answer []*SeedDev
	extendCurrSeed(seed, gg, read, true, false)
	//var newTStart, newQStart int32 = int32(seed.TargetStart) - 1, int32(seed.QueryStart) - 1
	if (seed.TargetStart <= 0 && seed.QueryStart > 0) && (len(gg.Nodes[seed.TargetId].Prev) > 0) {
		var depthSeeds []*SeedDev
		var prevLeft int
		for _, prev := range gg.Nodes[seed.TargetId].Prev {
			if len(prev.Dest.Seq) > 0 {
				var newTStart int32 = int32(len(prev.Dest.Seq)) - 1
				var newQStart int32 = int32(seed.QueryStart) - 1
				if prev.Dest.Seq[newTStart] == read.Seq[newQStart] {
					prevSeed := &SeedDev{TargetId: prev.Dest.Id, TargetStart: uint32(newTStart), QueryStart: uint32(newQStart), Length: 1, PosStrand: seed.PosStrand, Next: nil, Prev: nil}
					prevSeed.Next = seed
					depthSeeds = toTheLeft(prevSeed, gg, read)
					for prevLeft = 0; prevLeft < len(depthSeeds); prevLeft++ {
						answer = append(answer, depthSeeds[prevLeft])
					}
				}
			}
		}
	} else {
		answer = append(answer, seed)
	}
	return answer
}

func extendSeedTogether(seed *SeedDev, gg *SimpleGraph, read *fastq.Fastq) []*SeedDev {
	var answer []*SeedDev
	rightGraph := toTheRight(seed, gg, read)

	for rSeeds := 0; rSeeds < len(rightGraph); rSeeds++ {
		answer = append(answer, toTheLeft(rightGraph[rSeeds], gg, read)...)
	}
	return answer
}

func toTail(a *SeedDev) *SeedDev {
	if a.Next == nil {
		return a
	} else {
		return toTail(a.Next)
	}
}

func CompareSumLen(a *SeedDev, b *SeedDev) int {
	if sumLen(a) == sumLen(b) {
		return 0
	} else if sumLen(a) < sumLen(b) {
		return -1
	} else if sumLen(a) > sumLen(b) {
		return 1
	} else {
		log.Fatalf("Error: SeedDev len compare failed on:%d %d %d, %d %d %d\n", a.TargetId, a.TargetStart, a.Length, b.TargetId, b.TargetStart, b.Length)
		return 0
	}
}

func CompareBlastScore(a *SeedDev, b *SeedDev, read *fastq.Fastq) int {
	if BlastSeed(a, read) == BlastSeed(b, read) {
		return 0
	} else if BlastSeed(a, read) < BlastSeed(b, read) {
		return -1
	} else if BlastSeed(a, read) > BlastSeed(b, read) {
		return 1
	} else {
		log.Fatalf("Error: SeedDev len compare failed on:%d %d %d, %d %d %d\n", a.TargetId, a.TargetStart, a.Length, b.TargetId, b.TargetStart, b.Length)
		return 0
	}
}

func SortSeedExtended(seeds []*SeedDev) {
	sort.Slice(seeds, func(i, j int) bool { return CompareSumLen(seeds[i], seeds[j]) == 1 })
}

func SortBlastz(seeds []*SeedDev, read *fastq.Fastq) {
	sort.Slice(seeds, func(i, j int) bool { return CompareBlastScore(seeds[i], seeds[j], read) == 1 })
}

func BlastSeed(seed *SeedDev, read *fastq.Fastq) int64 {
	if seed.Next == nil {
		return scoreSeed(seed, read)
	} else {
		return scoreSeed(seed, read) + scoreSeed(seed.Next, read)
	}
}
/*
func IndexGenomeGraph(genome []*Node, seedLen int, seedStep int) map[uint64][]*SeedBed {
	if seedLen < 2 || seedLen > 32 {
		log.Fatalf("Error: seed length needs to be greater than 1 and less than 33.  Got: %d\n", seedLen)
	}
	answer := make(map[uint64][]*SeedBed)
	var seqCode uint64
	var nodeIdx, pos int
	for nodeIdx = 0; nodeIdx < len(genome); nodeIdx++ {
		for pos = 0; pos < len(genome[nodeIdx].Seq)-seedLen+1; pos += seedStep {
			if len(genome[nodeIdx].Seq)-pos >= seedLen {
				if dna.CountBaseInterval(genome[nodeIdx].Seq, dna.N, pos, pos+seedLen) == 0 {
					seqCode = dnaToNumber(genome[nodeIdx].Seq, pos, pos+seedLen)
					curr := SeedBed{Id: genome[nodeIdx].Id, Start: uint32(pos), End: uint32(pos + seedLen), Next: nil}
					answer[seqCode] = append(answer[seqCode], &curr)
					//if dna.CountBaseInterval(genome[nodeIdx].Seq, dna.N, pos, common.Min(pos+seedLen, len(genome[nodeIdx].Seq))) == 0 {
					//	seqCode = dnaToNumber(genome[nodeIdx].Seq, pos, pos+seedLen)
					//	curr := SeedBed{Id: genome[nodeIdx].Id, Start: uint32(pos), End: common.Min(pos+seedLen, len(genome[nodeIdx].Seq)), Next: nil}

					//	answer[seqCode] = append(answer[seqCode], &curr)

				}
			}

		}
	}
	return answer
}

func CheckSmallerHash(seedHash map[uint64][]*SeedBed, read *fastq.Fastq, seedLen int, stepSize int, posStrand bool, numMismatch int) []*SeedDev {
	var codedSeqs []uint64
	var hits []*SeedDev = make([]*SeedDev, 0)
	var i int
	for subSeqStart := 0; subSeqStart < len(read.Seq)-seedLen+1; subSeqStart += 2 {
		if dna.CountBaseInterval(read.Seq, dna.N, subSeqStart, subSeqStart+seedLen) == 0 {
			codedSeqs = maskSeedLen(read.Seq, subSeqStart, subSeqStart+seedLen, numMismatch)
			for i = 0; i < len(codedSeqs); i++ {
				currHits := seedHash[codedSeqs[i]]
				for _, value := range currHits {
					hits = append(hits, seedBedToSeedDev(value, uint32(subSeqStart), posStrand))
				}
			}
		}
	}
	log.Printf("Number of seeds found at 24 %d...\n", len(hits))
	return hits
}*/

func seedBedMask(a *SeedBed, currQPos uint32, posStrand bool, numMismatch int) *SeedDev {
	if a == nil {
		return nil
	} else {
		return &SeedDev{TargetId: a.Id, TargetStart: a.Start, QueryStart: currQPos, Length: a.End - a.Start - uint32(numMismatch), PosStrand: posStrand, Next: seedBedMask(a.Next, currQPos+a.End-a.Start-uint32(numMismatch), posStrand, numMismatch)}
	}
}
/*
func findSeedsInGraph(seedHash map[uint64][]*SeedBed, read *fastq.Fastq, seedLen int, stepSize int, posStrand bool, mask int) []*SeedDev {
	var codedSeq uint64 = 0
	var hits []*SeedDev = make([]*SeedDev, 0)
	for subSeqStart := 0; subSeqStart < len(read.Seq)-seedLen+1; subSeqStart += 2 {
		if dna.CountBaseInterval(read.Seq, dna.N, subSeqStart, subSeqStart+seedLen) == 0 {
			codedSeq = dnaToNumber(read.Seq, subSeqStart, subSeqStart+seedLen)
			currHits := seedHash[codedSeq]
			for _, value := range currHits {
				hits = append(hits, seedBedToSeedDev(value, uint32(subSeqStart), posStrand))
			}
		}
	}
	if len(hits) == 0 {
		log.Printf("No hits were found for seed length %d...\n", seedLen)
		log.Printf("Trying a smaller seed length...\n")
		hits = CheckSmallerHash(seedHash, read, seedLen, mask, posStrand, mask)
	}
	return hits
}*/

func seedBedToSeed(a *SeedBed, currQPos uint32, posStrand bool) *SeedDev {
	if a == nil {
		return nil
	} else {
		return &SeedDev{TargetId: a.Id, TargetStart: a.Start, QueryStart: currQPos, Length: a.End - a.Start, PosStrand: posStrand, Next: seedBedToSeedDev(a.Next, currQPos+a.End-a.Start, posStrand)}
	}
}

//TODO continue working on copying to head
/*
func extendSeedRight(seed *SeedDev, gg *SimpleGraph, read *fastq.Fastq) []*SeedDev {
	var newTEnd, newQEnd int32 = int32(seed.TargetStart + seed.Length), int32(seed.QueryStart + seed.Length)
	var graphGenomeHash []*SeedDev
	if int(newTEnd) == len(gg.Nodes[seed.TargetId].Seq) && int(newQEnd) < len(read.Seq) {
		var newTStart, newQStart int32
		//var numEdgeMatch int = 0
		for _, next := range gg.Nodes[seed.TargetId].Next {
			newTStart = 0
			newQStart = newQEnd
			if next.Dest.Seq[newTStart] == read.Seq[newQStart] {
				//numEdgeMatch++
				nextSeed := SeedDev{TargetId: next.Dest.Id, TargetStart: uint32(newTEnd), QueryStart: uint32(seed.QueryStart), Length: 1, PosStrand: seed.PosStrand, Next: nil, Prev: seed}
				if seed.Next != nil {
					seed.Next = &nextSeed
				} else {
					//copySeed := copySeed(seed)

					//copySeed = pointToHead(copySeed)
					//graphGenomeHash = append(graphGenomeHash, copySeed)
				}

				newTEnd, newQEnd = newTEnd+int32(nextSeed.Length), newQEnd+int32(nextSeed.Length)
				for ; int(newTEnd) < len(gg.Nodes[next.Dest.Id].Seq) && int(newQEnd) < len(read.Seq) && (gg.Nodes[seed.TargetId].Seq[newTEnd] == read.Seq[newQEnd]); newTEnd, newQEnd = newTEnd+1, newQEnd+1 {
					seed.Length++
				}
				extendSeedRight(&nextSeed, gg, read)
				//graphGenomeHash = append(graphGenomeHash, extendSeedRight(&nextSeed, gg, read)...)

			}
		}
	}
	graphGenomeHash = append(graphGenomeHash, seed)
	return graphGenomeHash
}
func extendSeedLeft(seed *SeedDev, gg *SimpleGraph, read *fastq.Fastq) []*SeedDev {
	var graphGenomeHash []*SeedDev
	if seed.TargetStart == 0 && seed.QueryStart > 0 {
		var newTStart, newQStart int32
		for _, prev := range gg.Nodes[seed.TargetId].Prev {
			newTStart = int32(len(prev.Dest.Seq)) - 1
			newQStart = int32(seed.QueryStart) - 1
			if prev.Dest.Seq[newTStart] == read.Seq[newQStart] {

				prevSeed := SeedDev{TargetId: prev.Dest.Id, TargetStart: uint32(newTStart), QueryStart: uint32(newQStart), Length: 1, PosStrand: seed.PosStrand, Next: nil, Prev: nil}
				prevSeed.Next = seed
				for newTStart, newQStart = newTStart-1, newQStart-1; newTStart >= 0 && newQStart >= 0 && (gg.Nodes[prevSeed.TargetId].Seq[newTStart] == read.Seq[newQStart]); newTStart, newQStart = newTStart-1, newQStart-1 {
					seed.TargetStart = uint32(newTStart)
					seed.QueryStart = uint32(newQStart)
					seed.Length++
				}
				graphGenomeHash = append(graphGenomeHash, extendSeedLeft(&prevSeed, gg, read)...)
			}
		}
	} else {
		graphGenomeHash = append(graphGenomeHash, seed)
	}

	return graphGenomeHash
}

func copyOfSeed(seed *SeedDev) *SeedDev {
	copyOfSeed := &SeedDev{TargetId: seed.TargetId, TargetStart: seed.TargetStart, QueryStart: seed.QueryStart, Length: seed.Length, PosStrand: seed.PosStrand, Next: nil, Prev: nil}
	return copyOfSeed
}

func pointToHead(seed *SeedDev) *SeedDev {
	if seed.Prev == nil {
		return seed
	} else {
		return pointToHead(seed.Prev)
	}
}*/
/*
func isNextSeedBetter(curr *Seed, currBestScore int64, perfectScore int64, queryLen int64, maxMatch int64, minMatch int64, leastSevereMismatch int64, leastSevereMatchMismatchChange int64) bool {
	seedLen := int64(FindTotalLengthSeed(curr, curr.Length))
	seeds := queryLen / (seedLen + 1)
	remainder := queryLen % (seedLen + 1)

	// seed by itself could be best
	if seedLen*maxMatch >= currBestScore &&
		perfectScore-((queryLen-seedLen)*minMatch) >= currBestScore {
		return true
		// seed along with whole seeds, but not the remainder
	} else if seedLen*seeds*maxMatch+seeds*leastSevereMismatch >= currBestScore &&
		perfectScore-remainder*minMatch+seeds*leastSevereMatchMismatchChange >= currBestScore {
		return true
		// seed along with whole seeds, as well as both remainders
	} else if seedLen*seeds*maxMatch+remainder*maxMatch+(seeds+1)*leastSevereMismatch >= currBestScore &&
		perfectScore+(seeds+1)*leastSevereMatchMismatchChange >= currBestScore {
		return true
	} else {
		return false
	}
}

func DevScoreSeed(seed *Seed, read *fastq.Fastq) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		score += HumanChimpTwoScoreMatrix[read.Seq[i]][read.Seq[i]]
		if seed.Next != nil {
			score += DevScoreSeed(seed.Next, read)
		}
	}
	return score
}

func SeedOverLap(a *Seed, b *Seed) bool {
	if a == nil || b == nil {
		return false
	}
	aLen := FindTotalLengthSeed(a, a.Length)
	bLen := FindTotalLengthSeed(b, b.Length)
	if common.MaxInt32(a.TargetStart, b.TargetStart) < common.MinInt32(a.TargetStart+aLen, b.TargetStart+bLen) {
		if common.MaxInt32(a.QueryStart, b.QueryStart) < common.MinInt32(a.QueryStart+aLen, b.QueryStart+bLen) {
			return true
		}
	}
	return false
}

func FindTotalLengthSeed(seed *Seed, length int32) int32 {
	if seed.Next != nil {
		length += seed.Next.Length
		FindTotalLengthSeed(seed.Next, length)
	}
	return length
}*/
