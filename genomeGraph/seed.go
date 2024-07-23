package genomeGraph

import (
	"container/heap"
	"log"
	"sort"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/numbers"
)

const basesPerInt int = 32

type Seed struct {
	TargetId    uint32
	TargetStart uint32
	QueryStart  uint32
	Length      uint32
	PosStrand   bool
	TotalLength uint32
	NextPart    *Seed
}

func extendToTheRight(node *Node, read *fastq.FastqBig, sk scoreKeeper, settings *GraphSettings, seedBuildHelper *SeedMemory, readStart int, nodeStart int, posStrand bool) []Seed {
	pq := make(PriorityQueue, 0, len(node.Next))
	heap.Init(&pq)

	seedBuildHelper.nodeOffset = nodeStart % basesPerInt
	seedBuildHelper.readOffset = 31 - ((readStart - seedBuildHelper.nodeOffset + 31) % 32)
	seedBuildHelper.rightMatches = 0
	var currNode Seed
	var nextParts []Seed
	var i, j int = 0, 0
	if posStrand {
		seedBuildHelper.rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, &read.Rainbow[seedBuildHelper.readOffset], readStart+seedBuildHelper.readOffset)
	} else {
		seedBuildHelper.rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, &read.RainbowRc[seedBuildHelper.readOffset], readStart+seedBuildHelper.readOffset)
	}
	// nothing aligned here
	if seedBuildHelper.rightMatches == 0 {
		return nil
	}

	settings.MaxMatch, settings.MinMatch, settings.LeastSevereMismatch, settings.LeastSevereMatchMismatchChange = MismatchStats(settings.ScoreMatrix)
	// we went all the way to end and there might be more
	if readStart+seedBuildHelper.rightMatches < len(read.Seq) && nodeStart+seedBuildHelper.rightMatches == node.SeqTwoBit.Len && len(node.Next) != 0 {
		for i = 0; i < len(node.Next); i++ {
			nextParts = extendToTheRight(node.Next[i].Dest, read, sk, settings, seedBuildHelper, readStart+seedBuildHelper.rightMatches, 0, posStrand)
			// if we aligned into the next node, make a seed for this node and point it to the next one
			for j = 0; j < len(nextParts); j++ {
				//currNode =
				seedBuildHelper.tempSeed = Seed{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(seedBuildHelper.rightMatches), PosStrand: posStrand, TotalLength: uint32(seedBuildHelper.rightMatches) + nextParts[j].TotalLength, NextPart: &nextParts[j]}

				heap.Push(&pq, &SeedHeap{Seed: &currNode})

			}
		}
	} else {
		currNode = Seed{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(seedBuildHelper.rightMatches), PosStrand: posStrand, TotalLength: uint32(seedBuildHelper.rightMatches), NextPart: nil}

		heap.Push(&pq, &SeedHeap{Seed: &currNode})

	}
	var answer []Seed
	// if the alignment did not go to another node, return the match for this node
	for pq.Len() > 0 {
		answer = append(answer, *heap.Pop(&pq).(*SeedHeap).Seed)
	}

	return answer
}

func extendToTheLeft(node *Node, read *fastq.FastqBig, seedBuildHelper *SeedMemory, currPart Seed) []Seed {
	var answers []Seed
	var i int
	var readBase dna.Base

	pq := make(PriorityQueue, 0)
	heap.Init(&pq)

	if currPart.QueryStart > 0 && currPart.TargetStart == 0 {
		for i = 0; i < len(node.Prev); i++ {
			if currPart.PosStrand {
				readBase = dnaTwoBit.GetBase(&read.Rainbow[0], uint(currPart.QueryStart)-1)
			} else {
				readBase = dnaTwoBit.GetBase(&read.RainbowRc[0], uint(currPart.QueryStart)-1)
			}
			if readBase == dnaTwoBit.GetBase(node.Prev[i].Dest.SeqTwoBit, uint(node.Prev[i].Dest.SeqTwoBit.Len)-1) {
				answers = extendToTheLeftHelper(node.Prev[i].Dest, read, seedBuildHelper, currPart, &pq)
				for _, part := range answers {
					heap.Push(&pq, &SeedHeap{Seed: &part})
				}
			}
		}
	} else {
		heap.Push(&pq, &SeedHeap{Seed: &currPart})
	}

	for pq.Len() > 0 {
		answers = append(answers, *heap.Pop(&pq).(*SeedHeap).Seed)
	}

	return answers
}

func extendToTheLeftHelper(node *Node, read *fastq.FastqBig, seedBuildHelper *SeedMemory, nextPart Seed, pq *PriorityQueue) []Seed {
	var nodePos int = node.SeqTwoBit.Len - 1
	var readPos int = int(nextPart.QueryStart) - 1
	seedBuildHelper.nodeOffset = nodePos % basesPerInt
	seedBuildHelper.readOffset = 31 - ((readPos - seedBuildHelper.nodeOffset + 31) % 32)
	seedBuildHelper.leftMatches = 0
	//var currPart Seed
	var answer []Seed
	var readBase dna.Base

	if nextPart.PosStrand {
		seedBuildHelper.leftMatches = numbers.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.Rainbow[seedBuildHelper.readOffset], readPos+seedBuildHelper.readOffset))
	} else {
		seedBuildHelper.leftMatches = numbers.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.RainbowRc[seedBuildHelper.readOffset], readPos+seedBuildHelper.readOffset))
	}

	if seedBuildHelper.leftMatches == 0 {
		log.Printf("Warning: should not have zero matches to the left\n")
		return nil // Indicate no match
	}
	seedBuildHelper.tempSeed = Seed{TargetId: node.Id, TargetStart: uint32(nodePos - (seedBuildHelper.leftMatches - 1)), QueryStart: uint32(readPos - (seedBuildHelper.leftMatches - 1)), Length: uint32(seedBuildHelper.leftMatches), PosStrand: nextPart.PosStrand, TotalLength: uint32(seedBuildHelper.leftMatches) + nextPart.TotalLength, NextPart: &nextPart}
	// we went all the way to end and there might be more
	if seedBuildHelper.tempSeed.QueryStart > 0 && seedBuildHelper.tempSeed.TargetStart == 0 {
		for _, prevNode := range node.Prev {
			if nextPart.PosStrand {
				readBase = dnaTwoBit.GetBase(&read.Rainbow[0], uint(seedBuildHelper.tempSeed.QueryStart)-1)
			} else {
				readBase = dnaTwoBit.GetBase(&read.RainbowRc[0], uint(seedBuildHelper.tempSeed.QueryStart)-1)
			}

			// Check if the previous node matches
			if readBase == dnaTwoBit.GetBase(prevNode.Dest.SeqTwoBit, uint(prevNode.Dest.SeqTwoBit.Len)-1) {
				answer = extendToTheLeftHelper(prevNode.Dest, read, seedBuildHelper, seedBuildHelper.tempSeed, pq)
				return append(answer, seedBuildHelper.tempSeed) // Append seedBuildHelper.tempSeed only once, after all recursive calls
			}
		}
	}
	// if the alignment did not go to another node, return the match for this node
	heap.Push(pq, &SeedHeap{Seed: &seedBuildHelper.tempSeed})
	return []Seed{seedBuildHelper.tempSeed}
}

// seedBuildHelper.nodeIdx, seedBuildHelper.nodePos int64 = 0, 0.
func seedMapMemPool(seedHash map[uint64][]uint64, nodes []Node, read *fastq.FastqBig, sk scoreKeeper, settings *GraphSettings, seedBuildHelper *SeedMemory) []Seed {
	const basesPerInt int64 = 32
	restartSeedHelper(seedBuildHelper)
	seedBuildHelper.keyShift = 64 - (uint(settings.TileSize) * 2)

	seedBuildHelper.Hits = seedBuildHelper.Hits[:0]
	seedBuildHelper.Worker = seedBuildHelper.Worker[:0]

	pq := make(PriorityQueue, 0, 20)
	heap.Init(&pq)

	for readStart := 0; readStart < len(read.Seq)-settings.TileSize+1; readStart++ {
		seedBuildHelper.keyIdx = (readStart + 31) / 32
		seedBuildHelper.keyOffset = 31 - ((readStart + 31) % 32)
		// do fwd strand
		seedBuildHelper.seqKey = read.Rainbow[seedBuildHelper.keyOffset].Seq[seedBuildHelper.keyIdx] >> seedBuildHelper.keyShift
		seedBuildHelper.currHits = seedHash[seedBuildHelper.seqKey]

		for _, seedBuildHelper.codedNodeCoord = range seedBuildHelper.currHits {
			seedBuildHelper.nodeIdx, seedBuildHelper.nodePos = numberToChromAndPos(seedBuildHelper.codedNodeCoord)
			seedBuildHelper.nodeOffset = int(seedBuildHelper.nodePos % basesPerInt)
			seedBuildHelper.readOffset = 31 - ((readStart - seedBuildHelper.nodeOffset + 31) % 32)
			seedBuildHelper.leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[seedBuildHelper.nodeIdx].SeqTwoBit, int(seedBuildHelper.nodePos), &read.Rainbow[seedBuildHelper.readOffset], readStart+seedBuildHelper.readOffset))
			seedBuildHelper.Worker = extendToTheRight(&nodes[seedBuildHelper.nodeIdx], read, sk, settings, seedBuildHelper, readStart-(seedBuildHelper.leftMatches-1), int(seedBuildHelper.nodePos)-(seedBuildHelper.leftMatches-1), true)
			for _, seedBuildHelper.tempSeed = range seedBuildHelper.Worker {
				for _, seed := range extendToTheLeft(&nodes[seedBuildHelper.nodeIdx], read, seedBuildHelper, seedBuildHelper.tempSeed) {
					heap.Push(&pq, &SeedHeap{Seed: &seed})
				}
				for pq.Len() > 0 {
					seedBuildHelper.Hits = append(seedBuildHelper.Hits, *heap.Pop(&pq).(*SeedHeap).Seed)
				}
			}
		}
		// do rev strand
		seedBuildHelper.seqKey = read.RainbowRc[seedBuildHelper.keyOffset].Seq[seedBuildHelper.keyIdx] >> seedBuildHelper.keyShift
		seedBuildHelper.currHits = seedHash[seedBuildHelper.seqKey]
		for _, seedBuildHelper.codedNodeCoord = range seedBuildHelper.currHits {
			seedBuildHelper.nodeIdx, seedBuildHelper.nodePos = numberToChromAndPos(seedBuildHelper.codedNodeCoord)
			seedBuildHelper.nodeOffset = int(seedBuildHelper.nodePos % basesPerInt)
			seedBuildHelper.readOffset = 31 - ((readStart - seedBuildHelper.nodeOffset + 31) % 32)

			seedBuildHelper.leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[seedBuildHelper.nodeIdx].SeqTwoBit, int(seedBuildHelper.nodePos), &read.RainbowRc[seedBuildHelper.readOffset], readStart+seedBuildHelper.readOffset))
			seedBuildHelper.Worker = extendToTheRight(&nodes[seedBuildHelper.nodeIdx], read, sk, settings, seedBuildHelper, readStart-(seedBuildHelper.leftMatches-1), int(seedBuildHelper.nodePos)-(seedBuildHelper.leftMatches-1), false)
			seedBuildHelper.Hits = append(seedBuildHelper.Hits, seedBuildHelper.Worker...)
		}
	}
	return seedBuildHelper.Hits
}

func getLastPart(a *Seed) *Seed {
	for ; a.NextPart != nil; a = a.NextPart {
	}
	return a
}

func SortSeedLen(seeds []Seed) {
	sort.Slice(seeds, func(i, j int) bool { return seeds[i].TotalLength > seeds[j].TotalLength })
}

func SortSeedDevByTotalLen(seeds []*Seed) {
	sort.Slice(seeds, func(i, j int) bool { return seeds[i].TotalLength > seeds[j].TotalLength })
}

func CompareBlastScore(a *Seed, b *Seed, read fastq.Fastq, scoreMatrix [][]int64) int {
	if BlastSeed(a, read, scoreMatrix) == BlastSeed(b, read, scoreMatrix) {
		return 0
	} else if BlastSeed(a, read, scoreMatrix) < BlastSeed(b, read, scoreMatrix) {
		return -1
	} else if BlastSeed(a, read, scoreMatrix) > BlastSeed(b, read, scoreMatrix) {
		return 1
	} else {
		log.Fatalf("Error: Seed len compare failed on:%d %d %d, %d %d %d\n", a.TargetId, a.TargetStart, a.Length, b.TargetId, b.TargetStart, b.Length)
		return 0
	}
}

func SortBlastz(seeds []*Seed, read fastq.Fastq, scoreMatrix [][]int64) {
	sort.Slice(seeds, func(i, j int) bool { return CompareBlastScore(seeds[i], seeds[j], read, scoreMatrix) == 1 })
}

func BlastSeed(seed *Seed, read fastq.Fastq, scoreMatrix [][]int64) int64 {
	if seed.NextPart == nil {
		return scoreSeed(seed, read, scoreMatrix)
	} else {
		return scoreSeed(seed, read, scoreMatrix) + scoreSeed(seed.NextPart, read, scoreMatrix)
	}
}
