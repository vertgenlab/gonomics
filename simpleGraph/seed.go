package simpleGraph

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

type Seed struct {
	Id  int64
	Start int64
	End int64
	//Score int64
}

func GraphTraversalHashMap(g *SimpleGraph, n *Node, idx int, seedLen int, ham5 map[uint64][]uint64) map[uint64][]uint64 {
	
	//answer := make(map[uint64][]uint64)
	var seqCode uint64
	if idx < len(n.Seq) - seedLen +1 {
		seqCode = dnaToNumber(n.Seq, idx, idx+seedLen)
		ham5[seqCode] = append(ham5[seqCode], chromAndPosToNumber(int(n.Id), idx))
		idx++
		ham5 = GraphTraversalHashMap(g, n, idx, seedLen, ham5)
	} else if (idx > len(n.Seq) - seedLen) && idx < len(n.Seq) {
		leftovers := make([]dna.Base, seedLen)
		copy(leftovers[0:len(n.Seq)-idx], n.Seq[idx:len(n.Seq)])
		
		for _, i := range n.Next {
			if len(i.Next.Seq) > seedLen {
				copy(leftovers[len(leftovers):seedLen], i.Next.Seq[0:seedLen-len(leftovers)])
				seqCode = dnaToNumber(leftovers, 0, seedLen)
				ham5[seqCode] = append(ham5[seqCode], chromAndPosToNumber(int(n.Id), idx))
				idx++
				ham5 = GraphTraversalHashMap(g, n, idx, seedLen, ham5)
			}
		}

	} else {
		for _, j := range n.Next {
			ham5 = GraphTraversalHashMap(g, j.Next, 0, seedLen, ham5)
		}
	}
	return ham5
}

func addSeed(existing []Seed, newId int64 , newStart int64 , newEnd int64) []Seed {
	//var answer
	var compareSeed int = -1
	for i:=0; i < len(existing) && compareSeed < 0; i++ {
		compareSeed = CompareSeed(existing[i], newId, newStart, newEnd)
		if  compareSeed < 0 {
			//do nothing

		//	existing = append(existing, Seed{Id: newId, Start: newStart, End: newEnd})
		} else if compareSeed == 0 {
			existing[i].Start = common.MinInt64(existing[i].Start, newStart)
			existing[i].End = common.MaxInt64(existing[i].End, newEnd)
			return existing
			//existing[i].Score++
			//need to check for more overlaps with i++
		} else if compareSeed > 0 {
			existing = append(existing, existing[len(existing)-1])
			copy(existing[i+1:], existing[i:])
			existing[i].Id = newId
			existing[i].Start = newStart
			existing[i].End = newEnd
			//existing[i].Score = 1
			return existing
		} else {
			log.Fatal("Some logic we did not catch when adding seed")
		}
	}
	existing = append(existing, Seed{Id: newId, Start: newStart, End: newEnd})
	return existing
}
//compare function
func CompareSeed(a Seed, newId int64 , newStart int64 , newEnd int64) int {
	if overlapSeed(a, newId, newStart, newEnd) {
		return 0
	} else if a.Id < newId || (a.Id == newId && a.Start < newStart) {
		return -1
	} else if a.Id > newId || (a.Id == newId && a.Start > newStart) {
		return 1
	} else {
		log.Fatal("Seed compare failed")
		return 0
	}
}

func overlapSeed(a Seed, newId int64, newStart int64, newEnd int64) bool {
	if a.Id == newId && common.MaxInt64(a.Start, newStart) < common.MinInt64(a.End, newEnd) {
		return true
	} else {
		return false
	}
}
