package genomeGraph

import (
	"log"
	"sort"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/numbers"
)

const (
	defaultMatrixSize int  = 2480
	leftTraversal     byte = 0
	rightTraversal    byte = 1
)

type SeedDev struct {
	TargetId    uint32
	TargetStart uint32
	QueryStart  uint32
	Length      uint32
	PosStrand   bool
	TotalLength uint32
	NextPart    *SeedDev
}

type seedHelper struct {
	currHits                                  []uint64
	codedNodeCoord                            uint64
	seqKey                                    uint64
	keyShift                                  uint
	keyIdx, keyOffset, readOffset, nodeOffset int
	nodeIdx, nodePos                          int64
	leftMatches                               int
	rightMatches                              int
	tempSeed                                  SeedDev
}

var HumanChimpTwoScoreMatrix = [][]int64{
	{90, -330, -236, -356, -208},
	{-330, 100, -318, -236, -196},
	{-236, -318, 100, -330, -196},
	{-356, -236, -330, 90, -208},
	{-208, -196, -196, -208, -202},
}

type ScoreMatrixHelper struct {
	Matrix                         [][]int64
	MaxMatch                       int64
	MinMatch                       int64
	LeastSevereMismatch            int64
	LeastSevereMatchMismatchChange int64
}

type memoryPool struct {
	Hits   []SeedDev
	Worker []SeedDev
}

type MatrixAln struct {
	m     [][]int64
	trace [][]byte
}

type dynamicScoreKeeper struct {
	i        int
	j        int
	routeIdx int
	currMax  int64
	route    []cigar.ByteCigar
}

type scoreKeeper struct {
	targetStart  int
	targetEnd    int
	queryStart   int
	queryEnd     int
	extension    int
	currScore    int64
	seedScore    int64
	perfectScore int64
	leftScore    int64
	rightScore   int64
	leftPath     []uint32
	rightPath    []uint32
	leftSeq      []dna.Base
	rightSeq     []dna.Base
	currSeq      []dna.Base
	tailSeed     SeedDev

	currSeed       SeedDev
	leftAlignment  []cigar.ByteCigar
	rightAlignment []cigar.ByteCigar
}

type dnaPool struct {
	Seq         []dna.Base
	Path        []uint32
	queryStart  int
	queryEnd    int
	targetStart int
	targetEnd   int
	currScore   int64
}

func NewDnaPool() sync.Pool {
	return sync.Pool{
		New: func() interface{} {
			dnaSeq := dnaPool{
				Seq:         make([]dna.Base, 0, 150),
				Path:        make([]uint32, 0, 10),
				queryStart:  0,
				targetStart: 0,
				targetEnd:   0,
				queryEnd:    0,
			}
			return &dnaSeq
		},
	}
}

func NewMemSeedPool() sync.Pool {
	return sync.Pool{
		New: func() interface{} {
			pool := memoryPool{
				Hits:   make([]SeedDev, 0, 10000),
				Worker: make([]SeedDev, 0, 10000),
			}
			return &pool
		},
	}
}

func createSeed(id uint32, readStart, nodeStart, length int, posStrand bool) *SeedDev {
	return &SeedDev{
		TargetId:    id,
		TargetStart: uint32(nodeStart),
		QueryStart:  uint32(readStart),
		Length:      uint32(length),
		PosStrand:   posStrand,
		TotalLength: uint32(length),
		NextPart:    nil,
	}
}

func extendToTheRight(node *Node, read fastq.FastqBig, readStart int, nodeStart int, posStrand bool) []*SeedDev {
	const basesPerInt int = 32
	var nodeOffset int = nodeStart % basesPerInt
	var readOffset int = 31 - ((readStart - nodeOffset + 31) % 32)
	var rightMatches int = 0
	var currNode *SeedDev = nil
	var answer, nextParts []*SeedDev
	var i, j int = 0, 0

	if posStrand {
		rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, &read.Rainbow[readOffset], readStart+readOffset)
	} else {
		rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, &read.RainbowRc[readOffset], readStart+readOffset)
	}

	// nothing aligned here
	if rightMatches == 0 {
		return nil
	}
	// we went all the way to end and there might be more
	if readStart+rightMatches < len(read.Seq) && nodeStart+rightMatches == node.SeqTwoBit.Len && len(node.Next) != 0 {
		for i = 0; i < len(node.Next); i++ {
			nextParts = extendToTheRight(node.Next[i].Dest, read, readStart+rightMatches, 0, posStrand)
			// if we aligned into the next node, make a seed for this node and point it to the next one
			for j = 0; j < len(nextParts); j++ {
				currNode = &SeedDev{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(rightMatches), PosStrand: posStrand, TotalLength: uint32(rightMatches) + nextParts[j].TotalLength, NextPart: nextParts[j]}
				answer = append(answer, currNode)
			}
		}
	}

	// if the alignment did not go to another node and the read sequence is not fully traversed, return the match for this node
	if len(answer) == 0 && readStart+rightMatches < len(read.Seq) {
		currNode = &SeedDev{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(rightMatches), PosStrand: posStrand, TotalLength: uint32(rightMatches), NextPart: nil}
		answer = []*SeedDev{currNode}
	}

	return answer
}

func extendToTheLeft(node *Node, read fastq.FastqBig, currPart *SeedDev) []*SeedDev {
	if currPart.QueryStart <= 0 || currPart.TargetStart != 0 {
		return []*SeedDev{currPart}
	}
	var readBase dna.Base
	var answer []*SeedDev
	for _, prev := range node.Prev {
		if currPart.PosStrand {
			readBase = dnaTwoBit.GetBase(&read.Rainbow[0], uint(currPart.QueryStart)-1)
		} else {
			readBase = dnaTwoBit.GetBase(&read.RainbowRc[0], uint(currPart.QueryStart)-1)
		}
		if readBase != dnaTwoBit.GetBase(prev.Dest.SeqTwoBit, uint(prev.Dest.SeqTwoBit.Len)-1) {
			continue
		}
		prevParts := extendToTheLeftHelper(prev.Dest, read, currPart)
		answer = append(answer, prevParts...)
	}
	return answer
}

func extendToTheLeftHelper(node *Node, read fastq.FastqBig, nextPart *SeedDev) []*SeedDev {
	const basesPerInt int = 32
	nodePos := node.SeqTwoBit.Len - 1
	readPos := int(nextPart.QueryStart) - 1
	nodeOffset := nodePos % basesPerInt
	readOffset := 31 - ((readPos - nodeOffset + 31) % 32)
	var leftMatches int
	if nextPart.PosStrand {
		leftMatches = numbers.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.Rainbow[readOffset], readPos+readOffset))
	} else {
		leftMatches = numbers.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, &read.RainbowRc[readOffset], readPos+readOffset))
	}
	if leftMatches == 0 {
		log.Fatal("Error: should not have zero matches to the left\n")
	}
	currPart := createSeed(uint32(node.Id), int(nodePos-(leftMatches-1)), int(readPos-(leftMatches-1)), int(uint32(leftMatches)+nextPart.TotalLength), nextPart.PosStrand)
	if currPart.QueryStart <= 0 || currPart.TargetStart != 0 {
		return []*SeedDev{currPart}
	}
	var answer []*SeedDev
	var readBase dna.Base
	for _, prev := range node.Prev {
		if nextPart.PosStrand {
			readBase = dnaTwoBit.GetBase(&read.Rainbow[0], uint(currPart.QueryStart)-1)
		} else {
			readBase = dnaTwoBit.GetBase(&read.RainbowRc[0], uint(currPart.QueryStart)-1)
		}
		if readBase != dnaTwoBit.GetBase(prev.Dest.SeqTwoBit, uint(prev.Dest.SeqTwoBit.Len)-1) {
			continue
		}
		prevParts := extendToTheLeftHelper(prev.Dest, read, currPart)
		answer = append(answer, prevParts...)
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

func findSeedsInSmallMapWithMemPool(seedHash map[uint64][]uint64, nodes []Node, read fastq.FastqBig, seedLen int, perfectScore int64, scoreMatrix [][]int64) []*SeedDev {
	const basesPerInt int64 = 32
	var currHits []uint64
	var codedNodeCoord uint64
	var leftMatches int = 0
	var keyIdx, keyOffset, readOffset, nodeOffset int = 0, 0, 0, 0
	var nodeIdx, nodePos int64 = 0, 0
	//var poolHead *SeedDev = *memoryPool
	var seqKey uint64
	var keyShift uint = 64 - (uint(seedLen) * 2)
	var tempSeeds, finalSeeds []*SeedDev
	var tempSeed *SeedDev

	for readStart := 0; readStart < len(read.Seq)-seedLen+1; readStart++ {
		keyIdx = (readStart + 31) / 32
		keyOffset = 31 - ((readStart + 31) % 32)

		// do fwd strand
		seqKey = read.Rainbow[keyOffset].Seq[keyIdx] >> keyShift
		currHits = seedHash[seqKey]
		/*if len(currHits) > 0 {
			log.Printf(" %d hits\n", len(currHits))
		}*/
		for _, codedNodeCoord = range currHits {
			nodeIdx, nodePos = numberToChromAndPos(codedNodeCoord)
			nodeOffset = int(nodePos % basesPerInt)
			readOffset = 31 - ((readStart - nodeOffset + 31) % 32)

			leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), &read.Rainbow[readOffset], readStart+readOffset))
			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.Rainbow[readOffset], readStart+readOffset)
			tempSeeds = extendToTheRight(&nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), true)
			//log.Printf("After extendToTheRight fwd:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				finalSeeds = append(finalSeeds, extendToTheLeft(&nodes[nodeIdx], read, tempSeed)...)
			}
			//log.Printf("After extendToTheLeft fwd:\n")
			//printSeedDev(finalSeeds)
			// TODO: Bring back speed optimizations once we are sure of correctness
			/*if leftMatches < 1 || rightMatches < seedLen {
				log.Fatalf("No matches found at seed location: %s %d, %d %d", dna.BasesToString(read.Seq), readStart, nodeIdx, nodePos)
			}*/

			/*if seedCouldBeBetter(int64(leftMatches+rightMatches-1), bestScore, perfectScore, int64(len(read.Seq)), 100, 90, -196, -296) {
				if poolHead != nil {
					currSeed = poolHead
					poolHead = poolHead.Next
				} else {
					currSeed = &SeedDev{}
				}
				currSeed.TargetId = uint32(nodeIdx)
				currSeed.TargetStart = uint32(int(nodePos) - leftMatches + 1)
				currSeed.QueryStart = uint32(readStart - leftMatches + 1)
				currSeed.Length = uint32(leftMatches + rightMatches - 1)
				currSeed.PosStrand = true
				currSeed.TotalLength = uint32(leftMatches + rightMatches - 1)
				currSeed.Next = hits
				hits = currSeed
				seedScore = scoreSeedFastqBig(currSeed, read, scoreMatrix)
				if seedScore > bestScore {
					bestScore = seedScore
				}
			}*/
		}

		// do rev strand
		seqKey = read.RainbowRc[keyOffset].Seq[keyIdx] >> keyShift
		currHits = seedHash[seqKey]
		/*if len(currHits) > 0 {
			log.Printf(" %d hits\n", len(currHits))
		}*/
		for _, codedNodeCoord = range currHits {
			nodeIdx, nodePos = numberToChromAndPos(codedNodeCoord)
			nodeOffset = int(nodePos % basesPerInt)
			readOffset = 31 - ((readStart - nodeOffset + 31) % 32)

			leftMatches = numbers.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), &read.RainbowRc[readOffset], readStart+readOffset))
			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.RainbowRc[readOffset], readStart+readOffset)
			tempSeeds = extendToTheRight(&nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), false)
			//log.Printf("After extendToTheRight rev:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				//log.Printf("tempSeed.QueryStart = %d\n", tempSeed.QueryStart)
				finalSeeds = append(finalSeeds, extendToTheLeft(&nodes[nodeIdx], read, tempSeed)...)
			}
			//log.Printf("After extendToTheLeft rev:\n")
			//printSeedDev(finalSeeds)
			// TODO: bring back speed optimizations
			/*if leftMatches < 1 || rightMatches < seedLen {
				log.Fatalf("No matches found at seed location: %s %d, %d %d", dna.BasesToString(read.SeqRc), readStart, nodeIdx, nodePos)
			}*/
			/*if seedCouldBeBetter(int64(leftMatches+rightMatches-1), bestScore, perfectScore, int64(len(read.SeqRc)), 100, 90, -196, -296) {
				if poolHead != nil {
					currSeed = poolHead
					poolHead = poolHead.Next
				} else {
					currSeed = &SeedDev{}
				}
				currSeed.TargetId = uint32(nodeIdx)
				currSeed.TargetStart = uint32(int(nodePos) - leftMatches + 1)
				currSeed.QueryStart = uint32(readStart - leftMatches + 1)
				currSeed.Length = uint32(leftMatches + rightMatches - 1)
				currSeed.PosStrand = false
				currSeed.TotalLength = uint32(leftMatches + rightMatches - 1)
				currSeed.Next = hits
				hits = currSeed
				seedScore = scoreSeedFastqBig(currSeed, read, scoreMatrix)
				if seedScore > bestScore {
					bestScore = seedScore
				}
			}*/
		}
	}

	/*var finalHits, nextSeed *SeedDev
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
	return hits*/
	//printSeedDev(finalSeeds)
	return finalSeeds
}
