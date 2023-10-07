package genomeGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"sort"
)

func extendCurrSeed(seed *SeedDev, gg *GenomeGraph, read fastq.Fastq, left bool, right bool) {
	var newTStart, newQStart, newTEnd, newQEnd int32 = int32(seed.TargetStart) - 1, int32(seed.QueryStart) - 1, int32(seed.TargetStart + seed.Length), int32(seed.QueryStart + seed.Length)
	//check to see if beginning is at index zero, if so do something like SeedDev.Prev
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

func toTheRight(seed *SeedDev, gg *GenomeGraph, read fastq.Fastq) []*SeedDev {
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
					nextSeed := &SeedDev{TargetId: next.Dest.Id, TargetStart: uint32(newTStart), QueryStart: uint32(newQStart), Length: 1, PosStrand: seed.PosStrand, NextPart: nil}
					edgeSeeds = toTheRight(nextSeed, gg, read)
					for e = 0; e < len(edgeSeeds); e++ {
						currSeed := &SeedDev{TargetId: seed.TargetId, TargetStart: seed.TargetStart, QueryStart: seed.QueryStart, Length: seed.Length, PosStrand: seed.PosStrand, NextPart: edgeSeeds[e]}
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

func toTheLeft(seed *SeedDev, gg *GenomeGraph, read fastq.Fastq) []*SeedDev {
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
					prevSeed := &SeedDev{TargetId: prev.Dest.Id, TargetStart: uint32(newTStart), QueryStart: uint32(newQStart), Length: 1, PosStrand: seed.PosStrand, NextPart: seed}
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

func extendSeedTogether(seed *SeedDev, gg *GenomeGraph, read fastq.Fastq) []*SeedDev {
	var answer []*SeedDev
	rightGraph := toTheRight(seed, gg, read)

	for rSeeds := 0; rSeeds < len(rightGraph); rSeeds++ {
		answer = append(answer, toTheLeft(rightGraph[rSeeds], gg, read)...)
	}
	return answer
}

func getLastPart(a *SeedDev) *SeedDev {
	for ; a.NextPart != nil; a = a.NextPart {
	}
	return a
}

func CompareBlastScore(a *SeedDev, b *SeedDev, read fastq.Fastq, scoreMatrix [][]int64) int {
	scoreA := BlastSeed(a, read, scoreMatrix)
	scoreB := BlastSeed(b, read, scoreMatrix)

	switch {
	case scoreA < scoreB:
		return -1
	case scoreA > scoreB:
		return 1
	default:
		return 0
	}
}

var blastSeedMemo map[*SeedDev]int64 = make(map[*SeedDev]int64)

func BlastSeed(seed *SeedDev, read fastq.Fastq, scoreMatrix [][]int64) int64 {
	// Check if the score is already computed for the seed
	if val, exists := blastSeedMemo[seed]; exists {
		return val
	}

	// Compute the score and save it for future use
	var score int64
	if seed.NextPart == nil {
		score = scoreSeed(seed, read, scoreMatrix)
	} else {
		score = scoreSeed(seed, read, scoreMatrix) + scoreSeed(seed.NextPart, read, scoreMatrix)
	}
	blastSeedMemo[seed] = score
	return score
}

func SortBlastz(seeds []*SeedDev, read fastq.Fastq, scoreMatrix [][]int64) {
	sort.Slice(seeds, func(i, j int) bool { return CompareBlastScore(seeds[i], seeds[j], read, scoreMatrix) == 1 })
}
