package simpleGraph

import (
	"github.com/vertgenlab/gonomics/dna"
	//	"github.com/vertgenlab/gonomics/fastq"
	"log"
	"sort"
)

type SeedDev struct {
	TargetId    uint32
	TargetStart uint32
	QueryStart  uint32
	Length      uint32
	PosStrand   bool
	TotalLength uint32
	NextPart    *SeedDev
	Next        *SeedDev
}

func IndexGenomeIntoMap(genome []*Node, seedLen int, seedStep int) map[uint64][]uint64 {
	if seedLen < 2 || seedLen > 32 {
		log.Fatalf("Error: seed length needs to be greater than 1 and less than 33.  Got: %d\n", seedLen)
	}
	answer := make(map[uint64][]uint64)
	var seqCode, locationCode uint64
	var nodeIdx, pos int
	for nodeIdx = 0; nodeIdx < len(genome); nodeIdx++ {
		for pos = 0; pos < len(genome[nodeIdx].Seq)-seedLen+1; pos += seedStep {
			if dna.CountBaseInterval(genome[nodeIdx].Seq, dna.N, pos, pos+seedLen) == 0 {
				seqCode = dnaToNumber(genome[nodeIdx].Seq, pos, pos+seedLen)
				answer[seqCode] = append(answer[seqCode], ChromAndPosToNumber(nodeIdx, pos))
			}
		}
		for ; pos < len(genome[nodeIdx].Seq); pos += seedStep {
			locationCode = ChromAndPosToNumber(nodeIdx, pos)
			for edgeIdx := 0; edgeIdx < len(genome[nodeIdx].Next); edgeIdx++ {
				indexGenomeIntoMapHelper(genome[nodeIdx].Seq[pos:], genome[nodeIdx].Next[edgeIdx].Dest, locationCode, seedLen, answer)
			}
		}
	}
	return answer
}

func indexGenomeIntoMapHelper(prevSeq []dna.Base, currNode *Node, locationCode uint64, seedLen int, seedMap map[uint64][]uint64) {
	if len(prevSeq)+len(currNode.Seq) >= seedLen {
		currSeq := append(prevSeq, currNode.Seq[0:(seedLen-len(prevSeq))]...)
		if dna.CountBaseInterval(currSeq, dna.N, 0, seedLen) == 0 {
			seqCode := dnaToNumber(currSeq, 0, seedLen)
			seedMap[seqCode] = append(seedMap[seqCode], locationCode)
		}
	} else {
		for edgeIdx := 0; edgeIdx < len(currNode.Next); edgeIdx++ {
			indexGenomeIntoMapHelper(append(prevSeq, currNode.Seq...), currNode.Next[edgeIdx].Dest, locationCode, seedLen, seedMap)
		}
	}
}

func indexGenomeIntoSlice(genome []*Node, seedLen int, seedStep int) [][]uint64 {
	if seedLen < 2 || seedLen > 16 {
		log.Fatalf("Error: seed length needs to be greater than 1 and less than 17.  Got: %d\n", seedLen)
	}
	var sliceSize int
	sliceSize = 1 << uint(seedLen*2)
	answer := make([][]uint64, sliceSize)
	var seqCode, locationCode uint64
	var nodeIdx, pos int
	//var seqScratch []dna.Base = make([]dna.Base, seedLen)
	for nodeIdx = 0; nodeIdx < len(genome); nodeIdx++ {
		for pos = 0; pos < len(genome[nodeIdx].Seq)-seedLen+1; pos += seedStep {
			if dna.CountBaseInterval(genome[nodeIdx].Seq, dna.N, pos, pos+seedLen) == 0 {
				seqCode = dnaToNumber(genome[nodeIdx].Seq, pos, pos+seedLen)
				answer[seqCode] = append(answer[seqCode], ChromAndPosToNumber(nodeIdx, pos))
			}
		}
		for ; pos < len(genome[nodeIdx].Seq); pos += seedStep {
			locationCode = ChromAndPosToNumber(nodeIdx, pos)
			for edgeIdx := 0; edgeIdx < len(genome[nodeIdx].Next); edgeIdx++ {
				indexGenomeIntoSliceHelper(genome[nodeIdx].Seq[pos:], genome[nodeIdx].Next[edgeIdx].Dest, locationCode, seedLen, answer)
			}
		}
	}
	return answer
}

func indexGenomeIntoSliceHelper(prevSeq []dna.Base, currNode *Node, locationCode uint64, seedLen int, seedSlice [][]uint64) {
	if len(prevSeq)+len(currNode.Seq) >= seedLen {
		currSeq := append(prevSeq, currNode.Seq[0:(seedLen-len(prevSeq))]...)
		if dna.CountBaseInterval(currSeq, dna.N, 0, seedLen) == 0 {
			seqCode := dnaToNumber(currSeq, 0, seedLen)
			seedSlice[seqCode] = append(seedSlice[seqCode], locationCode)
		}
	} else {
		for edgeIdx := 0; edgeIdx < len(currNode.Next); edgeIdx++ {
			indexGenomeIntoSliceHelper(append(prevSeq, currNode.Seq...), currNode.Next[edgeIdx].Dest, locationCode, seedLen, seedSlice)
		}
	}
}

// TODO: this does not take into account breaking up seeds by gaps instead of mismatches
// similar calculations could also be used as the parameters to a banded alignment
func seedCouldBeBetter(seedLen int64, currBestScore int64, perfectScore int64, queryLen int64, maxMatch int64, minMatch int64, leastSevereMismatch int64, leastSevereMatchMismatchChange int64) bool {
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

func numberOfSeeds(head *SeedDev) int {
	var answer int = 0
	for ; head != nil; head = head.Next {
		answer++
	}
	return answer
}

func SortSeedDevListByTotalLen(headPtr **SeedDev) {
	var head *SeedDev = *headPtr
	if (head != nil) && (head.Next != nil) {
		var halfway *SeedDev = nil
		halfway = SplitSeedDevListInHalf(head)
		SortSeedDevListByTotalLen(&head)
		SortSeedDevListByTotalLen(&halfway)
		*headPtr = MergeSortedSeedDevs(head, halfway)
	}
}

func MergeSortedSeedDevs(a *SeedDev, b *SeedDev) *SeedDev {
	var dummyHead *SeedDev = &SeedDev{}
	var curr *SeedDev = dummyHead

	for (a != nil) && (b != nil) {
		if a.TotalLength >= b.TotalLength {
			curr.Next = a
			curr = curr.Next
			a = a.Next
		} else {
			curr.Next = b
			curr = curr.Next
			b = b.Next
		}
	}
	if a == nil && b != nil {
		curr.Next = b
	}
	if b == nil && a != nil {
		curr.Next = a
	}
	return dummyHead.Next
}

func SplitSeedDevListInHalf(head *SeedDev) *SeedDev {
	var fast, slow, secondHalf *SeedDev = nil, nil, nil
	slow = head
	fast = head.Next

	for fast != nil {
		fast = fast.Next
		if fast != nil {
			slow = slow.Next
			fast = fast.Next
		}
	}
	secondHalf = slow.Next
	slow.Next = nil
	return secondHalf
}

func SortSeedDevByTotalLen(seeds []*SeedDev) {
	sort.Slice(seeds, func(i, j int) bool { return seeds[i].TotalLength > seeds[j].TotalLength })
}
