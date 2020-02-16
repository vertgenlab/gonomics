package simpleGraph

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
	"sort"
)

type SeedBed struct {
	Id    uint32
	Start uint32
	End   uint32
	Next  *SeedBed
}

type SeedDev struct {
	TargetId    uint32
	TargetStart uint32
	QueryStart  uint32
	Length      uint32
	PosStrand   bool
	Next        *SeedDev
	Prev        *SeedDev
}

func extendSeedsDev(seeds []*SeedDev, gg *SimpleGraph, read *fastq.Fastq) {
	for i := 0; i < len(seeds); i++ {
		extendSeedDev(seeds[i], gg, read)
	}
}

// TODO: only works for fasta style graphs, no following edges
func extendSeedDev(seed *SeedDev, gg *SimpleGraph, read *fastq.Fastq) {
	var newTStart, newQStart, newTEnd, newQEnd int32 = int32(seed.TargetStart) - 1, int32(seed.QueryStart) - 1, int32(seed.TargetStart + seed.Length), int32(seed.QueryStart + seed.Length)
	//check to see if begining is at index zero, if so do something like SeedDev.Prev
	//if newStart < 0
	for ; newTStart >= 0 && newQStart >= 0 && (gg.Nodes[seed.TargetId].Seq[newTStart] == read.Seq[newQStart]); newTStart, newQStart = newTStart-1, newQStart-1 {
		seed.TargetStart = uint32(newTStart)
		seed.QueryStart = uint32(newQStart)
		seed.Length++
	}
	for ; int(newTEnd) < len(gg.Nodes[seed.TargetId].Seq) && int(newQEnd) < len(read.Seq) && (gg.Nodes[seed.TargetId].Seq[newTEnd] == read.Seq[newQEnd]); newTEnd, newQEnd = newTEnd+1, newQEnd+1 {
		seed.Length++
	}
}

func printSeedDev(a []*SeedDev) {
	for i, _ := range a {
		log.Printf("%d\t%d\t%d\t%d\t%t\n", a[i].TargetId, a[i].TargetStart, a[i].QueryStart, a[i].Length, a[i].PosStrand)
	}
}

func sumLen(a *SeedDev) uint32 {
	if a.Next == nil {
		return a.Length
	} else {
		return a.Length + sumLen(a.Next)
	}
}

func SortIndex(tiles [][]*SeedBed) {
	for i := 0; i < len(tiles); i++ {
		SortSeedBedByPos(tiles[i])
	}
}

func ReduceIndexCapacity(tiles [][]*SeedBed) {
	for i := 0; i < len(tiles); i++ {
		tmp := make([]*SeedBed, len(tiles[i]))
		copy(tmp, tiles[i])
		tiles[i] = tmp
	}
}

func indexGenomeHelper(n *Node, seedLen int) []*SeedBed {
	if len(n.Seq) >= seedLen {
		if dna.CountBaseInterval(n.Seq, dna.N, 0, seedLen) != 0 {
			return []*SeedBed{}
		} else {
			curr := SeedBed{Id: n.Id, Start: 0, End: uint32(seedLen)}
			return []*SeedBed{&curr}
		}
	} else {
		if dna.CountBaseInterval(n.Seq, dna.N, 0, len(n.Seq)) != 0 || len(n.Next) == 0 {
			return []*SeedBed{}
		} else {
			answer := make([]*SeedBed, 0)
			for i := 0; i < len(n.Next); i++ {
				endings := indexGenomeHelper(n.Next[i].Dest, seedLen-len(n.Seq))
				for j := 0; j < len(endings); j++ {
					curr := SeedBed{Id: n.Id, Start: 0, End: uint32(len(n.Seq)), Next: endings[j]}
					answer = append(answer, &curr)
				}
			}
			return answer
		}
	}
}

func IndexGenomeIntoMap(genome []*Node, seedLen int, seedStep int) map[uint64][]*SeedBed {
	if seedLen < 2 || seedLen > 32 {
		log.Fatalf("Error: seed length needs to be greater than 1 and less than 33.  Got: %d\n", seedLen)
	}
	answer := make(map[uint64][]*SeedBed)
	var seqCode uint64
	var nodeIdx, pos int
	for nodeIdx = 0; nodeIdx < len(genome); nodeIdx++ {
		for pos = 0; pos < len(genome[nodeIdx].Seq)-seedLen+1; pos += seedStep {
			if dna.CountBaseInterval(genome[nodeIdx].Seq, dna.N, pos, pos+seedLen) == 0 {
				seqCode = dnaToNumber(genome[nodeIdx].Seq, pos, pos+seedLen)
				curr := SeedBed{Id: genome[nodeIdx].Id, Start: uint32(pos), End: uint32(pos + seedLen), Next: nil}
				answer[seqCode] = append(answer[seqCode], &curr)
			}
		}
	}
	return answer
}

func IndexGenomeIntoSlice(genome []*Node, seedLen int, seedStep int) [][]*SeedBed {
	if seedLen < 2 || seedLen > 16 {
		log.Fatalf("Error: seed length needs to be greater than 1 and less than 16.  Got: %d\n", seedLen)
	}
	var sliceSize int
	sliceSize = 1 << uint(seedLen*2)
	//fmt.Printf("Slice size is :%d, seedLen=%d\n", sliceSize, seedLen)
	answer := make([][]*SeedBed, sliceSize)
	var seqCode uint64
	var nodeIdx, pos int
	for nodeIdx = 0; nodeIdx < len(genome); nodeIdx++ {
		for pos = 0; pos < len(genome[nodeIdx].Seq)-seedLen+1; pos += seedStep {
			if dna.CountBaseInterval(genome[nodeIdx].Seq, dna.N, pos, pos+seedLen) == 0 {
				seqCode = dnaToNumber(genome[nodeIdx].Seq, pos, pos+seedLen)
				curr := SeedBed{Id: genome[nodeIdx].Id, Start: uint32(pos), End: uint32(pos + seedLen), Next: nil}
				answer[seqCode] = append(answer[seqCode], &curr)
			}
		}
	}
	return answer
}

// TODO: this does not take into account breaking up seeds by gaps instead of mismatches
// similar calculations could also be used as the parameters to a banded alignment
func seedCouldBeBetter(curr *SeedDev, currBestScore int64, perfectScore int64, queryLen int64, maxMatch int64, minMatch int64, leastSevereMismatch int64, leastSevereMatchMismatchChange int64) bool {
	seedLen := int64(curr.Length)
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

/*
func addSeedDev(existing []*SeedDev, curr *SeedDev) []*SeedDev {
	var compareSeed int = -1
	for i := 0; i < len(existing) && compareSeed < 0; i++ {
		if canMerge(existing[i], curr) {
			existing[i].TargetStart = common.MinInt64(existing[i].TargetStart, curr.TargetStart)
			existing[i].TargetEnd = common.MaxInt64(existing[i].TargetEnd, curr.TargetEnd)
			existing[i].QueryStart = common.MinInt64(existing[i].QueryStart, curr.QueryStart)
			existing[i].QueryEnd = common.MaxInt64(existing[i].QueryEnd, curr.QueryEnd)
			return existing
		}
		compareSeed = CompareSeedDev(existing[i], curr)
		if compareSeed <= 0 {
			//do nothing
		} else if compareSeed > 0 {
			existing = append(existing, existing[len(existing)-1])
			copy(existing[i+1:], existing[i:])
			existing[i] = curr
			return existing
		} else {
			log.Fatal("Some logic we did not catch when adding seed")
		}
	}
	existing = append(existing, curr)
	return existing
}*/

func canMerge(a *SeedDev, b *SeedDev) bool {
	if a.TargetId == b.TargetId &&
		common.MaxUint32(a.TargetStart, b.TargetStart) <= common.MinUint32(a.TargetStart+a.Length, b.TargetStart+b.Length) &&
		common.MaxUint32(a.QueryStart, b.QueryStart) <= common.MinUint32(a.QueryStart+a.Length, b.QueryStart+b.Length) &&
		int(b.TargetStart)-int(a.TargetStart) == int(b.QueryStart)-int(a.QueryStart) {
		return true
	} else {
		return false
	}
}

func SortSeedBedByPos(seeds []*SeedBed) {
	sort.Slice(seeds, func(i, j int) bool { return CompareSeedBed(seeds[i], seeds[j]) == -1 })
}

func CompareSeedBed(a *SeedBed, b *SeedBed) int {
	if a.Id < b.Id ||
		(a.Id == b.Id && a.Start < b.Start) ||
		(a.Id == b.Id && a.Start == b.Start && a.End < b.End) {
		return -1
	} else if a.Id > b.Id ||
		(a.Id == b.Id && a.Start > b.Start) ||
		(a.Id == b.Id && a.Start == b.Start && a.End > b.End) {
		return 1
	} else if a.Id == b.Id && a.Start == b.Start && a.End == b.End {
		return 0
	} else {
		log.Fatalf("Error: SeedBed compare failed on:%d %d %d, %d %d %d\n", a.Id, a.Start, a.End, b.Id, b.Start, b.End)
		return 0
	}
}

func SortSeedDevByLen(seeds []*SeedDev) {
	sort.Slice(seeds, func(i, j int) bool { return CompareLenSeedDev(seeds[i], seeds[j]) == 1 })
}

func CompareLenSeedDev(a *SeedDev, b *SeedDev) int {
	if a.Length == b.Length {
		return 0
	} else if a.Length < b.Length {
		return -1
	} else if a.Length > b.Length {
		return 1
	} else {
		log.Fatalf("Error: SeedDev len compare failed on:%d %d %d, %d %d %d\n", a.TargetId, a.TargetStart, a.Length, b.TargetId, b.TargetStart, b.Length)
		return 0
	}
}

func CompareSeedDev(a *SeedDev, b *SeedDev) int {
	if a.TargetId == b.TargetId && a.TargetStart == b.TargetStart && a.Length == b.Length {
		return 0
	} else if a.TargetId < b.TargetId ||
		(a.TargetId == b.TargetId && a.TargetStart < b.TargetStart) ||
		(a.TargetId == b.TargetId && a.TargetStart == b.TargetStart && a.Length < b.Length) {
		return -1
	} else if a.TargetId > b.TargetId ||
		(a.TargetId == b.TargetId && a.TargetStart > b.TargetStart) ||
		(a.TargetId == b.TargetId && a.TargetStart == b.TargetStart && a.Length > b.Length) {
		return 1
	} else {
		log.Fatalf("Error: SeedDev compare failed on:%d %d %d, %d %d %d\n", a.TargetId, a.TargetStart, a.Length, b.TargetId, b.TargetStart, b.Length)
		return 0
	}
}
