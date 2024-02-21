package genomeGraph

import (
	"log"
	"sort"

	"github.com/vertgenlab/gonomics/fastq"
)

func extendCurrSeed(seed *SeedDev, gg *GenomeGraph, read fastq.Fastq, extendLeft, extendRight bool) {
	nodeSeq := gg.Nodes[seed.TargetId].Seq
	readSeq := read.Seq

	if extendLeft {
		for newTStart, newQStart := int(seed.TargetStart)-1, int(seed.QueryStart)-1; newTStart >= 0 && newQStart >= 0 && nodeSeq[newTStart] == readSeq[newQStart]; newTStart, newQStart = newTStart-1, newQStart-1 {
			seed.TargetStart = uint32(newTStart)
			seed.QueryStart = uint32(newQStart)
			seed.Length++
		}
	}

	if extendRight {
		for newTEnd, newQEnd := int(seed.TargetStart+seed.Length), int(seed.QueryStart+seed.Length); newTEnd < len(nodeSeq) && newQEnd < len(readSeq) && nodeSeq[newTEnd] == readSeq[newQEnd]; newTEnd, newQEnd = newTEnd+1, newQEnd+1 {
			seed.Length++
		}
	}
}

func toTheRight(seed *SeedDev, gg *GenomeGraph, read fastq.Fastq) []*SeedDev {
	extendCurrSeed(seed, gg, read, false, true)
	return extendToDirection(seed, gg, read, true)
}

func toTheLeft(seed *SeedDev, gg *GenomeGraph, read fastq.Fastq) []*SeedDev {
	extendCurrSeed(seed, gg, read, true, false)
	return extendToDirection(seed, gg, read, false)
}

func extendToDirection(seed *SeedDev, gg *GenomeGraph, read fastq.Fastq, toRight bool) []*SeedDev {
	var answer []*SeedDev
	if toRight {
		// Existing right extension logic...
	} else {
		// Extending to the left
		targetStart := int(seed.TargetStart)
		queryStart := int(seed.QueryStart)

		if targetStart <= 0 && queryStart > 0 && len(gg.Nodes[seed.TargetId].Prev) > 0 {
			for _, edge := range gg.Nodes[seed.TargetId].Prev {
				prevNodeEnd := len(edge.Dest.Seq) - 1
				prevNodeQueryStart := queryStart - 1

				if edge.Dest.Seq[prevNodeEnd] == read.Seq[prevNodeQueryStart] {
					prevSeed := &SeedDev{
						TargetId:    edge.Dest.Id,
						TargetStart: uint32(prevNodeEnd),
						QueryStart:  uint32(prevNodeQueryStart),
						Length:      1,
						PosStrand:   seed.PosStrand,
						NextPart:    seed,
					}
					depthSeeds := toTheLeft(prevSeed, gg, read)
					for _, depthSeed := range depthSeeds {
						answer = append(answer, depthSeed)
					}
				}
			}
		} else {
			answer = append(answer, seed)
		}
	}
	return answer
}

func extendSeedTogether(seed *SeedDev, gg *GenomeGraph, read fastq.Fastq) []*SeedDev {
	rightSeeds := toTheRight(seed, gg, read)
	var answer []*SeedDev

	for _, rightSeed := range rightSeeds {
		leftSeeds := toTheLeft(rightSeed, gg, read)
		answer = append(answer, leftSeeds...)
	}

	return answer
}

func getLastPart(a *SeedDev) *SeedDev {
	for ; a.NextPart != nil; a = a.NextPart {
	}
	return a
}

func CompareBlastScore(a *SeedDev, b *SeedDev, read fastq.Fastq, scoreMatrix [][]int64) int {
	if BlastSeed(a, read, scoreMatrix) == BlastSeed(b, read, scoreMatrix) {
		return 0
	} else if BlastSeed(a, read, scoreMatrix) < BlastSeed(b, read, scoreMatrix) {
		return -1
	} else if BlastSeed(a, read, scoreMatrix) > BlastSeed(b, read, scoreMatrix) {
		return 1
	} else {
		log.Fatalf("Error: SeedDev len compare failed on:%d %d %d, %d %d %d\n", a.TargetId, a.TargetStart, a.Length, b.TargetId, b.TargetStart, b.Length)
		return 0
	}
}

func SortSeedsByScore(seeds []*SeedDev, read fastq.Fastq, scoreMatrix [][]int64) {
	sort.Slice(seeds, func(i, j int) bool {
		return BlastSeed(seeds[i], read, scoreMatrix) > BlastSeed(seeds[j], read, scoreMatrix)
	})
}

func SortBlastz(seeds []*SeedDev, read fastq.Fastq, scoreMatrix [][]int64) {
	sort.Slice(seeds, func(i, j int) bool { return CompareBlastScore(seeds[i], seeds[j], read, scoreMatrix) == 1 })
}

func BlastSeed(seed *SeedDev, read fastq.Fastq, scoreMatrix [][]int64) int64 {
	if seed.NextPart == nil {
		return scoreSeed(seed, read, scoreMatrix)
	} else {
		return scoreSeed(seed, read, scoreMatrix) + scoreSeed(seed.NextPart, read, scoreMatrix)
	}
}
