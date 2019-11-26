package simpleGraph

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
)

/*
type TournamentNode struct {
	CurrValue  *SeedDev
	WaitingValues []*SeedBed
	QueryStart uint32
	Left   *tournamentNode
	Right  *tournamentNode
}

func buildTree(leaves []*TournamentTree) *TournamentTree {
	if len(leaves) == 0 {
		log.Errorf("Error: tree building should not get down to a node of zero leaves\n")
		return nil
	} else if len(leaves) == 1 {
		var first *SeedBed = nil
		first, leaves[0].WaitingValues = leaves[0].WaitingValues[0], leaves[0].WaitingValues[1:]
		leaves[0].CurrValue = &SeedDev{TargetId: first.Id, TargetStart: first.Start, QueryStart: leaves[0], Length: first.Length
		return &leaves[0]
	} else {
		midpoint := len(leaves) / 2
		return &TournamentTree{currValue: nil, waitingValues: []*SeedBed{}, queryStart: 0, left:buildTree(leaves[0:midpoint]), right: buildTree(leaves[midpoint:len(leaves)])}
	}
}

*/

func seedBedToSeedDev(a *SeedBed, currQPos uint32, posStrand bool) *SeedDev {
	if a == nil {
		return nil
	} else {
		return &SeedDev{TargetId: a.Id, TargetStart: a.Start, QueryStart: currQPos, Length: a.End - a.Start, PosStrand: posStrand, Next: seedBedToSeedDev(a.Next, currQPos+a.End-a.Start, posStrand)}
	}
}

func mergeIntoSeed(aTail *SeedDev, b *SeedBed, currQPos uint32, posStrand bool) bool {
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
}

func getTails(a []*SeedDev) []*SeedDev {
	tails := make([]*SeedDev, len(a))
	var tail *SeedDev = nil
	for i := 0; i < len(a); i++ {
		for tail = a[0]; tail.Next != nil; tail = tail.Next {
		}
		tails[i] = tail
	}
	return tails
}

func mergeSeedLists(lastPosition []*SeedDev, currPosition []*SeedBed, currQPos uint32, posStrand bool) ([]*SeedDev, []*SeedDev) {
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
}

// need to handle neg strand
func findSeedsFast(seedHash [][]*SeedBed, read *fastq.Fastq, seedLen int) []*SeedDev {
	var codedSeq uint64 = 0
	var prevHits []*SeedDev = make([]*SeedDev, 0)
	var allHits []*SeedDev = make([]*SeedDev, 0)
	for subSeqStart := 0; subSeqStart < len(read.Seq)-seedLen+1; subSeqStart++ {
		codedSeq = dnaToNumber(read.Seq, subSeqStart, subSeqStart+seedLen)
		currHits := seedHash[codedSeq]
		noMerge, merged := mergeSeedLists(prevHits, currHits, uint32(subSeqStart), true)
		allHits = append(allHits, noMerge...)
		prevHits = merged
	}
	allHits = append(allHits, prevHits...)
	return allHits
}
