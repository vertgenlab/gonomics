package simpleGraph

import (
	//"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
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

func IndexGenomeGraph(genome []*Node, seedLen int, seedStep int) map[uint64][]*SeedBed {
	if seedLen < 2 || seedLen > 32 {
		log.Fatalf("Error: seed length needs to be greater than 1 and less than 33.  Got: %d\n", seedLen)
	}
	answer := make(map[uint64][]*SeedBed)
	var seqCode uint64
	var nodeIdx, pos int
	for nodeIdx = 0; nodeIdx < len(genome); {
		for pos = 0; pos < len(genome[nodeIdx].Seq)-seedLen+1; pos += 1 {
			if dna.CountBaseInterval(genome[nodeIdx].Seq, dna.N, pos, pos+seedLen) == 0 {
				seqCode = dnaToNumber(genome[nodeIdx].Seq, pos, pos+seedLen)
				curr := SeedBed{Id: genome[nodeIdx].Id, Start: uint32(pos), End: uint32(pos + seedLen), Next: nil}
				answer[seqCode] = append(answer[seqCode], &curr)
			}
		}
	}
	return answer
}

func genomeGraphDictionary(seeds []*SeedDev, gg *SimpleGraph, read *fastq.Fastq) []*SeedDev {
	var graphGenomeHash []*SeedDev

	for i := 0; i < len(seeds); i++ {
		graphGenomeHash = append(graphGenomeHash, extendSeedTogether(seeds[i], gg, read)...)
	}
	return graphGenomeHash
}

func findSeedsInGraph(seedHash map[uint64][]*SeedBed, read *fastq.Fastq, seedLen int, stepSize int, posStrand bool) []*SeedDev {
	var codedSeq uint64 = 0
	var hits []*SeedDev = make([]*SeedDev, 0)
	var curr *SeedDev
	for subSeqStart := 0; subSeqStart < len(read.Seq)-seedLen+1; subSeqStart++ {
		if dna.CountBaseInterval(read.Seq, dna.N, subSeqStart, subSeqStart+seedLen) == 0 {
			codedSeq = dnaToNumber(read.Seq, subSeqStart, subSeqStart+seedLen)
			currHits := seedHash[codedSeq]
			for _, value := range currHits {
				curr = seedBedToSeed(value, uint32(subSeqStart), posStrand)
				hits = append(hits, curr)
			}
		}
	}
	return hits
}

func seedBedToSeed(a *SeedBed, currQPos uint32, posStrand bool) *SeedDev {
	if a == nil {
		return nil
	} else {
		return &SeedDev{TargetId: a.Id, TargetStart: a.Start, QueryStart: currQPos, Length: a.End - a.Start, PosStrand: posStrand, Next: seedBedToSeedDev(a.Next, currQPos+a.End-a.Start, posStrand)}
	}
}

func extendSeedTogether(seed *SeedDev, gg *SimpleGraph, read *fastq.Fastq) []*SeedDev {
	var graphGenomeHash []*SeedDev
	extendSeedDev(seed, gg, read)
	rightGraph := extendSeedRight(seed, gg, read)
	for toTheLeft := 0; toTheLeft < len(rightGraph); toTheLeft++ {
		graphGenomeHash = append(graphGenomeHash, extendSeedLeft(rightGraph[toTheLeft], gg, read)...)
	}
	return graphGenomeHash
}
//TODO continue working on copying to head
func extendSeedRight(seed *SeedDev, gg *SimpleGraph, read *fastq.Fastq) []*SeedDev {
	var graphGenomeHash []*SeedDev
	var newTEnd, newQEnd int32 = int32(seed.TargetStart + seed.Length), int32(seed.QueryStart + seed.Length)
	if int(newTEnd) == len(gg.Nodes[seed.TargetId].Seq) && int(newQEnd) < len(read.Seq) {
		var newTStart, newQStart int32
		for _, next := range gg.Nodes[seed.TargetId].Next {
			newTEnd = 0
			newQStart = newQEnd
			if next.Dest.Seq[newTStart] == read.Seq[newQStart] {
				nextSeed := SeedDev{TargetId: next.Dest.Id, TargetStart: uint32(newTEnd), QueryStart: uint32(seed.QueryStart), Length: 1, PosStrand: seed.PosStrand, Next: nil, Prev: nil}
				nextSeed.Prev = copySeedToHead(seed)
				//graphGenomeHash = append(graphGenomeHash, seed)
				newTEnd, newQEnd = newTEnd+int32(nextSeed.Length), newQEnd+int32(nextSeed.Length)
				for ; int(newTEnd) < len(gg.Nodes[next.Dest.Id].Seq) && int(newQEnd) < len(read.Seq) && (gg.Nodes[seed.TargetId].Seq[newTEnd] == read.Seq[newQEnd]); newTEnd, newQEnd = newTEnd+1, newQEnd+1 {
					seed.Length++
				}
				graphGenomeHash = append(graphGenomeHash, extendSeedRight(&nextSeed, gg, read)...)
			}
		}
	}
	return graphGenomeHash
}
func extendSeedLeft(seed *SeedDev, gg *SimpleGraph, read *fastq.Fastq) []*SeedDev {
	var graphGenomeHash []*SeedDev
	if seed.TargetStart == 0 && seed.QueryStart > 0 {
		//make a new copy of yourself if prev already exists?
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
	}

	return graphGenomeHash
}

func copySeedToHead(seed *SeedDev) *SeedDev {
	copyOfPrev := &SeedDev{TargetId: seed.Prev.TargetId, TargetStart: seed.Prev.TargetStart, QueryStart: seed.Prev.QueryStart, Length: seed.Prev.Length, PosStrand: seed.Prev.PosStrand, Next: nil, Prev: nil}
	//var copyOfCurr SeedDev =
	if seed.Prev != nil {
		copyOfPrev = copySeedToHead(seed.Prev)
	}
	return copyOfPrev
}
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
