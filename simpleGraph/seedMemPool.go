package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaTwoBit"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"log"
	"sort"
	"sync"
	//"os"
	//"runtime"
	//"runtime/pprof"
)

type pathSeed struct {
	QStart int
	QPos   bool
	TStart int
	TEnd   int
	TNodes []uint32
	Len    int
}
type SeedPool struct {
	//answer  []pathSeed
	working []*pathSeed
}

func RoutinePathSeed(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, input <-chan FastqGsw, output chan<- GirafGsw, wg *sync.WaitGroup) {
	matrix := NewSwMatrix(10000)
	scorekeeper := scoreKeeper{}
	//seedPool := NewSeedPool()
	var seedPool = sync.Pool{
		New: func() interface{} {
			return make([]*SeedDev, 0, 1000)
		},
	}
	dynamicKeeper := dynamicScoreKeeper{}
	for read := range input {
		output <- wrapPathSeedGsw(gg, read, seedHash, seedLen, stepSize, &matrix, scoreMatrix, scorekeeper, dynamicKeeper, &seedPool)
	}
	wg.Done()
}
func wrapPathSeedGsw(gg *SimpleGraph, fq FastqGsw, seedHash map[uint64][]uint64, seedLen int, stepSize int, matrix *MatrixAln, scoreMatrix [][]int64, sk scoreKeeper, dynamicScore dynamicScoreKeeper, seedPool *sync.Pool) GirafGsw {
	var mappedPair GirafGsw = GirafGsw{
		ReadOne: pathSeedDevGraphSmithWaterman(gg, fq.ReadOne, seedHash, seedLen, stepSize, matrix, scoreMatrix, sk, dynamicScore, seedPool),
		ReadTwo: pathSeedDevGraphSmithWaterman(gg, fq.ReadTwo, seedHash, seedLen, stepSize, matrix, scoreMatrix, sk, dynamicScore, seedPool),
	}
	//setGirafFlags(&mappedPair)
	return mappedPair
}
func sortPathSeed(seeds []*pathSeed) {
	sort.Slice(seeds, func(i, j int) bool { return lenPathSeedCompare(seeds[i], seeds[j]) == 1 })
}

func lenPathSeedCompare(a *pathSeed, b *pathSeed) int {
	if a.Len == b.Len {
		return 0
	} else if a.Len < b.Len {
		return -1
	} else if a.Len > b.Len {
		return 1
	} else {
		log.Fatalf("Error: SeedDev total length compare failed...\n")
		return 0
	}
}

func pathSeedDevGraphSmithWaterman(gg *SimpleGraph, read fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, matrix *MatrixAln, scoreMatrix [][]int64, sk scoreKeeper, dynamicScore dynamicScoreKeeper, seedPool *sync.Pool) giraf.Giraf {
	var currBest giraf.Giraf = giraf.Giraf{
		QName:     read.Name,
		QStart:    0,
		QEnd:      0,
		PosStrand: true,
		Path:      &giraf.Path{},
		Cigar:     nil,
		AlnScore:  0,
		MapQ:      255,
		Seq:       read.Seq,
		Qual:      read.Qual,
		Notes:     []giraf.Note{giraf.Note{Tag: "XO", Type: 'Z', Value: "~"}},
	}
	resetScoreKeeper(sk)
	sk.perfectScore = perfectMatchBig(&read, scoreMatrix)
	sk.extension = int(sk.perfectScore/600) + len(read.Seq)
	//seedsHash := seedPool.Get().([]*pathSeed)
	//tempSeeds := extendSeeds.Get().([]*pathSeed)
	//var tmpSeeds []*SeedDev = extendPool.Get().([]*SeedDev)
	seeds := GetSeedPathHits(seedHash, gg.Nodes, &read, seedLen, sk.perfectScore, scoreMatrix, seedPool)
	//tempSeeds = tempSeeds[:0]
	//extendSeeds.Put(tempSeeds)
	sortPathSeed(seeds)
	//var tailSeed *SeedDev
	//var seedScore int64
	var currSeq []dna.Base = make([]dna.Base, len(read.Seq))
	var currSeed *pathSeed

	//leftSeq,

	var dnaPool = sync.Pool{
		New: func() interface{} {
			return make([]dna.Base, sk.extension)
		},
	}
	for i := 0; i < len(seeds) && seedCouldBeBetter(int64(seeds[i].Len), int64(currBest.AlnScore), sk.perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		currSeed = seeds[i]
		//tailSeed = getLastPart(currSeed)
		if currSeed.QPos {
			currSeq = read.Seq
		} else {
			currSeq = read.SeqRc
		}
		sk.seedScore = scoreSeedSeq(currSeq, uint32(currSeed.QStart), uint32(currSeed.QStart+currSeed.Len), scoreMatrix)
		if int(currSeed.Len) == len(currSeq) {
			sk.currScore = sk.seedScore
			sk.minTarget = currSeed.TStart
			sk.maxTarget = currSeed.TStart + currSeed.Len
			sk.minQuery = currSeed.QStart
			//sk.maxQuery = currSeed.Len - 1
		} else {
			sk.leftAlignment, sk.leftScore, sk.minTarget, sk.minQuery, sk.leftPath = LeftAlignTraversal(gg.Nodes[currSeed.TNodes[0]], sk.leftSeq, currSeed.TStart, sk.leftPath, sk.extension-currSeed.Len, currSeq[:currSeed.QStart], scoreMatrix, matrix, dynamicScore, &dnaPool)
			sk.rightAlignment, sk.rightScore, sk.maxTarget, sk.maxQuery, sk.rightPath = RightAlignTraversal(gg.Nodes[currSeed.TNodes[len(currSeed.TNodes)-1]], sk.rightSeq, currSeed.TEnd, sk.rightPath, sk.extension-currSeed.Len, currSeq[currSeed.QStart+currSeed.Len:], scoreMatrix, matrix, dynamicScore, &dnaPool)
			sk.currScore = sk.leftScore + sk.seedScore + sk.rightScore
		}
		if sk.currScore > int64(currBest.AlnScore) {
			currBest.QStart = sk.minQuery
			currBest.QEnd = currSeed.QStart + currSeed.Len + sk.maxQuery - sk.minQuery - 1
			currBest.PosStrand = currSeed.QPos
			ReversePath(sk.leftPath)
			currBest.Path = setPath(currBest.Path, sk.minTarget, CatPaths(CatPaths(sk.leftPath, currSeed.TNodes), sk.rightPath), sk.maxTarget)
			currBest.Cigar = SoftClipBases(sk.minQuery, len(currSeq), cigar.CatByteCigar(cigar.AddCigarByte(sk.leftAlignment, cigar.ByteCigar{RunLen: uint16(currSeed.Len), Op: 'M'}), sk.rightAlignment))
			currBest.AlnScore = int(sk.currScore)
			currBest.Seq = currSeq
			if &gg.Nodes[currBest.Path.Nodes[0]].Info != nil {
				currBest.Notes[0].Value = fmt.Sprintf("%s=%d", gg.Nodes[currBest.Path.Nodes[0]].Name, gg.Nodes[currBest.Path.Nodes[0]].Info.Start)
				currBest.Notes = append(currBest.Notes, infoToNotes(gg.Nodes, currBest.Path.Nodes))
			} else {
				currBest.Notes[0].Value = fmt.Sprintf("%s=%d", gg.Nodes[currBest.Path.Nodes[0]].Name, 1)
			}
			//if currBest.QEnd > len(currBest.Seq) {
			///	log.Fatalf("Error in calculating start and end query (minQuery=%d, seedStart=%d, seedLen=%d, maxQuery=%d)\n", sk.minQuery, currSeed.QStart, currSeed.Len, sk.maxQuery)
			//}
		}
	}
	//seeds = seeds[:0]
	//seedPool.Put(seeds)
	if !currBest.PosStrand {
		fastq.ReverseQualUint8Record(currBest.Qual)
	}
	return currBest
}

func GetSeedPathHits(seedHash map[uint64][]uint64, nodes []*Node, read *fastq.FastqBig, seedLen int, perfectScore int64, scoreMatrix [][]int64, seedPool *sync.Pool) []*pathSeed {
	const basesPerInt int64 = 32
	var currHits []uint64
	var codedNodeCoord uint64
	var leftMatches int = 0
	var keyIdx, keyOffset, readOffset, nodeOffset int = 0, 0, 0, 0
	var nodeIdx, nodePos int64 = 0, 0
	var finalSeeds []*pathSeed
	//var poolHead *SeedDev = *memoryPool
	var seqKey uint64
	var keyShift uint = 64 - (uint(seedLen) * 2)

	//var tempSeeds []*pathSeed //seedPool.Get().([]*SeedDev)
	var tempSeed *pathSeed

	pool := seedPool.Get().(*SeedPool)
	pool.working = pool.working[:0]
	for readStart := 0; readStart < len(read.Seq)-seedLen+1; readStart++ {
		keyIdx = (readStart + 31) / 32
		keyOffset = 31 - ((readStart + 31) % 32)

		// do fwd strand
		seqKey = read.Rainbow[keyOffset].Seq[keyIdx] >> keyShift
		currHits = seedHash[seqKey]

		for _, codedNodeCoord = range currHits {
			nodeIdx, nodePos = numberToChromAndPos(codedNodeCoord)
			nodeOffset = int(nodePos % basesPerInt)
			readOffset = 31 - ((readStart - nodeOffset + 31) % 32)

			leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.Rainbow[readOffset], readStart+readOffset))
			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.Rainbow[readOffset], readStart+readOffset)

			pool.working = ExtendRightPath(nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), true, seedPool)
			//log.Printf("After extendToTheRightDev fwd:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range pool.working {

				finalSeeds = append(finalSeeds, ExtendLeftPath(nodes[nodeIdx], read, tempSeed)...)
			}

		}
		// do rev strand
		seqKey = read.RainbowRc[keyOffset].Seq[keyIdx] >> keyShift
		currHits = seedHash[seqKey]
		/*if len(currHits) > 0 {
			log.Printf(" %d hits\n", len(currHits))
		}*/
		pool.working = pool.working[:0]
		for _, codedNodeCoord = range currHits {
			nodeIdx, nodePos = numberToChromAndPos(codedNodeCoord)
			nodeOffset = int(nodePos % basesPerInt)
			readOffset = 31 - ((readStart - nodeOffset + 31) % 32)

			leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.RainbowRc[readOffset], readStart+readOffset))

			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.RainbowRc[readOffset], readStart+readOffset)
			pool.working = ExtendRightPath(nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), false, seedPool)
			//log.Printf("After extendToTheRightDev rev:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range pool.working {
				//log.Printf("tempSeed.QueryStart = %d\n", tempSeed.QueryStart)
				finalSeeds = append(finalSeeds, ExtendLeftPath(nodes[nodeIdx], read, tempSeed)...)
			}

		}
	}
	seedPool.Put(pool)
	return finalSeeds
}

func ExtendRightPath(node *Node, read *fastq.FastqBig, readStart int, nodeStart int, posStrand bool, seedPool *sync.Pool) []*pathSeed {

	var answer []*pathSeed
	const basesPerInt int = 32
	var nodeOffset int = nodeStart % basesPerInt
	var readOffset int = 31 - ((readStart - nodeOffset + 31) % 32)
	var rightMatches int = 0
	var currNode *pathSeed
	currNode.TStart = nodeStart
	currNode.QStart = readStart
	currNode.QPos = posStrand
	currNode.TNodes = []uint32{node.Id}

	//var nextParts []pathSeed
	var i, j int = 0, 0

	if posStrand {
		rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, read.Rainbow[readOffset], readStart+readOffset)
	} else {
		rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, read.RainbowRc[readOffset], readStart+readOffset)
	}
	// nothing aligned here
	if rightMatches == 0 {
		return nil
	}
	// we went all the way to end and there might be more
	if readStart+rightMatches < len(read.Seq) && nodeStart+rightMatches == node.SeqTwoBit.Len && len(node.Next) != 0 {
		for i = 0; i < len(node.Next); i++ {
			//nextParts := extendSeeds.Get().([]*SeedDev)
			pool := seedPool.Get().(*SeedPool)
			pool.working = pool.working[:0]
			pool.working = ExtendRightPath(node.Next[i].Dest, read, readStart+rightMatches, 0, posStrand, seedPool)
			//nextParts =
			// if we aligned into the next node, make a seed for this node and point it to the next one
			for j = 0; j < len(pool.working); j++ {

				currNode.TNodes = append(currNode.TNodes, pool.working[j].TNodes...)
				currNode.TEnd = pool.working[j].TEnd
				currNode.Len = pool.working[j].Len
				answer = append(answer, currNode)
			}
			seedPool.Put(pool)
		}
	}
	// if the alignment did not go to another node, return the match for this node
	if len(answer) == 0 {
		currNode.Len += rightMatches
		currNode.TEnd = currNode.TStart + rightMatches
		answer = []*pathSeed{currNode}
	}
	return answer
}

func ExtendLeftPath(node *Node, read *fastq.FastqBig, currPart *pathSeed) []*pathSeed {
	var prevParts, answer []*pathSeed
	var i int
	var readBase dna.Base
	if currPart.QStart > 0 && currPart.TStart == 0 {
		for i = 0; i < len(node.Prev); i++ {
			if currPart.QPos {
				readBase = dnaTwoBit.GetBase(read.Rainbow[0], uint(currPart.QStart)-1)
			} else {
				readBase = dnaTwoBit.GetBase(read.RainbowRc[0], uint(currPart.QStart)-1)
			}
			if readBase == dnaTwoBit.GetBase(node.Prev[i].Dest.SeqTwoBit, uint(node.Prev[i].Dest.SeqTwoBit.Len)-1) {
				prevParts = leftHelperDev(node.Prev[i].Dest, read, currPart, prevParts)
				answer = append(answer, prevParts...)
			}
		}
	}
	if len(answer) == 0 {
		return []*pathSeed{currPart}
	} else {
		return answer
	}
}

func leftHelperDev(node *Node, read *fastq.FastqBig, nextPart *pathSeed, answer []*pathSeed) []*pathSeed {
	answer = answer[:0]
	const basesPerInt int = 32
	var nodePos int = node.SeqTwoBit.Len - 1
	var readPos int = int(nextPart.QStart) - 1
	var nodeOffset int = nodePos % basesPerInt
	var readOffset int = 31 - ((readPos - nodeOffset + 31) % 32)
	var leftMatches int = 0
	var currPart *pathSeed
	var prevParts []*pathSeed
	var i int
	var readBase dna.Base
	if nextPart.QPos {
		leftMatches = common.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, read.Rainbow[readOffset], readPos+readOffset))
	} else {
		leftMatches = common.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, read.RainbowRc[readOffset], readPos+readOffset))
	}
	if leftMatches == 0 {
		log.Fatal("Error: should not have zero matches to the left\n")
	}
	currPart = &pathSeed{QStart: readPos - (leftMatches - 1), QPos: nextPart.QPos, TStart: nodePos - (leftMatches - 1), TNodes: nextPart.TNodes, Len: leftMatches + nextPart.Len}
	// we went all the way to end and there might be more
	if currPart.QStart > 0 && currPart.TStart == 0 {
		for i = 0; i < len(node.Prev); i++ {
			if nextPart.QPos {
				readBase = dnaTwoBit.GetBase(read.Rainbow[0], uint(currPart.QStart)-1)
			} else {
				readBase = dnaTwoBit.GetBase(read.RainbowRc[0], uint(currPart.QStart)-1)
			}
			if readBase == dnaTwoBit.GetBase(node.Prev[i].Dest.SeqTwoBit, uint(node.Prev[i].Dest.SeqTwoBit.Len)-1) {
				prevParts = leftHelperDev(node.Prev[i].Dest, read, currPart, prevParts)
				answer = append(answer, prevParts...)
			}
		}
	}
	// if the alignment did not go to another node, return the match for this node
	if len(answer) == 0 {
		answer = []*pathSeed{currPart}
	}
	return answer
}

func extendToTheRightDev(node *Node, read *fastq.FastqBig, readStart int, nodeStart int, posStrand bool, answer []*SeedDev) []*SeedDev {
	const basesPerInt int = 32
	answer = answer[:0]
	var nodeOffset int = nodeStart % basesPerInt
	var readOffset int = 31 - ((readStart - nodeOffset + 31) % 32)
	var rightMatches int = 0
	var currNode SeedDev
	var nextParts []*SeedDev
	var i, j int = 0, 0

	if posStrand {
		rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, read.Rainbow[readOffset], readStart+readOffset)
	} else {
		rightMatches = dnaTwoBit.CountRightMatches(node.SeqTwoBit, nodeStart, read.RainbowRc[readOffset], readStart+readOffset)
	}
	// nothing aligned here
	if rightMatches == 0 {
		return nil
	}
	// we went all the way to end and there might be more
	if readStart+rightMatches < len(read.Seq) && nodeStart+rightMatches == node.SeqTwoBit.Len && len(node.Next) != 0 {
		for i = 0; i < len(node.Next); i++ {
			//nextParts := extendSeeds.Get().([]*SeedDev)
			nextParts = extendToTheRightDev(node.Next[i].Dest, read, readStart+rightMatches, 0, posStrand, nextParts)
			// if we aligned into the next node, make a seed for this node and point it to the next one
			for j = 0; j < len(nextParts); j++ {
				currNode = SeedDev{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(rightMatches), PosStrand: posStrand, TotalLength: uint32(rightMatches) + nextParts[j].TotalLength, NextPart: nextParts[j], Next: nil}
				answer = append(answer, &currNode)
			}
			//nextParts = nextParts[:0]
			//extendSeeds.Put(nextParts)
		}
	}
	// if the alignment did not go to another node, return the match for this node
	if len(answer) == 0 {
		currNode = SeedDev{TargetId: node.Id, TargetStart: uint32(nodeStart), QueryStart: uint32(readStart), Length: uint32(rightMatches), PosStrand: posStrand, TotalLength: uint32(rightMatches), NextPart: nil, Next: nil}
		answer = []*SeedDev{&currNode}
	}

	return answer
}

func extendToTheLeftDev(node *Node, read *fastq.FastqBig, currPart *SeedDev) []*SeedDev {
	var answer, prevParts []*SeedDev
	var i int
	var readBase dna.Base

	if currPart.QueryStart > 0 && currPart.TargetStart == 0 {
		for i = 0; i < len(node.Prev); i++ {
			if currPart.PosStrand {
				readBase = dnaTwoBit.GetBase(read.Rainbow[0], uint(currPart.QueryStart)-1)
			} else {
				readBase = dnaTwoBit.GetBase(read.RainbowRc[0], uint(currPart.QueryStart)-1)
			}
			if readBase == dnaTwoBit.GetBase(node.Prev[i].Dest.SeqTwoBit, uint(node.Prev[i].Dest.SeqTwoBit.Len)-1) {
				prevParts = extendToTheLeftHelperDev(node.Prev[i].Dest, read, currPart)
				answer = append(answer, prevParts...)
			}
		}
	}

	if len(answer) == 0 {
		return []*SeedDev{currPart}
	} else {
		return answer
	}
}

func extendToTheLeftHelperDev(node *Node, read *fastq.FastqBig, nextPart *SeedDev) []*SeedDev {
	const basesPerInt int = 32
	var nodePos int = node.SeqTwoBit.Len - 1
	var readPos int = int(nextPart.QueryStart) - 1
	var nodeOffset int = nodePos % basesPerInt
	var readOffset int = 31 - ((readPos - nodeOffset + 31) % 32)
	var leftMatches int = 0
	var currPart *SeedDev = nil
	var prevParts, answer []*SeedDev
	var i int
	var readBase dna.Base

	if nextPart.PosStrand {
		leftMatches = common.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, read.Rainbow[readOffset], readPos+readOffset))
	} else {
		leftMatches = common.Min(readPos+1, dnaTwoBit.CountLeftMatches(node.SeqTwoBit, nodePos, read.RainbowRc[readOffset], readPos+readOffset))
	}

	if leftMatches == 0 {
		log.Fatal("Error: should not have zero matches to the left\n")
	}

	currPart = &SeedDev{TargetId: node.Id, TargetStart: uint32(nodePos - (leftMatches - 1)), QueryStart: uint32(readPos - (leftMatches - 1)), Length: uint32(leftMatches), PosStrand: nextPart.PosStrand, TotalLength: uint32(leftMatches) + nextPart.TotalLength, NextPart: nextPart, Next: nil}

	// we went all the way to end and there might be more
	if currPart.QueryStart > 0 && currPart.TargetStart == 0 {
		for i = 0; i < len(node.Prev); i++ {
			if nextPart.PosStrand {
				readBase = dnaTwoBit.GetBase(read.Rainbow[0], uint(currPart.QueryStart)-1)
			} else {
				readBase = dnaTwoBit.GetBase(read.RainbowRc[0], uint(currPart.QueryStart)-1)
			}
			if readBase == dnaTwoBit.GetBase(node.Prev[i].Dest.SeqTwoBit, uint(node.Prev[i].Dest.SeqTwoBit.Len)-1) {
				prevParts = extendToTheLeftHelperDev(node.Prev[i].Dest, read, currPart)
				answer = append(answer, prevParts...)
			}
		}
	}

	// if the alignment did not go to another node, return the match for this node
	if len(answer) == 0 {
		answer = []*SeedDev{currPart}
	}
	return answer
}

func seedMapMemPool(seedHash map[uint64][]uint64, nodes []*Node, read *fastq.FastqBig, seedLen int, perfectScore int64, scoreMatrix [][]int64, finalSeeds []*SeedDev, tempSeeds []*SeedDev) []*SeedDev {
	const basesPerInt int64 = 32
	var currHits []uint64
	var codedNodeCoord uint64
	var leftMatches int = 0
	var keyIdx, keyOffset, readOffset, nodeOffset int = 0, 0, 0, 0
	var nodeIdx, nodePos int64 = 0, 0
	//var poolHead *SeedDev = *memoryPool
	var seqKey uint64
	var keyShift uint = 64 - (uint(seedLen) * 2)
	//var tempSeeds []*SeedDev
	//tempSeeds := seedPool.Get().([]*SeedDev)
	var tempSeed *SeedDev
	for readStart := 0; readStart < len(read.Seq)-seedLen+1; readStart++ {
		keyIdx = (readStart + 31) / 32
		keyOffset = 31 - ((readStart + 31) % 32)

		// do fwd strand
		seqKey = read.Rainbow[keyOffset].Seq[keyIdx] >> keyShift
		currHits = seedHash[seqKey]

		for _, codedNodeCoord = range currHits {
			nodeIdx, nodePos = numberToChromAndPos(codedNodeCoord)
			nodeOffset = int(nodePos % basesPerInt)
			readOffset = 31 - ((readStart - nodeOffset + 31) % 32)

			leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.Rainbow[readOffset], readStart+readOffset))
			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.Rainbow[readOffset], readStart+readOffset)
			tempSeeds = extendToTheRightDev(nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), true, tempSeeds)
			//log.Printf("After extendToTheRightDev fwd:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				finalSeeds = append(finalSeeds, extendToTheLeftDev(nodes[nodeIdx], read, tempSeed)...)
			}
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

			leftMatches = common.Min(readStart+1, dnaTwoBit.CountLeftMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.RainbowRc[readOffset], readStart+readOffset))
			//rightMatches = dnaTwoBit.CountRightMatches(nodes[nodeIdx].SeqTwoBit, int(nodePos), read.RainbowRc[readOffset], readStart+readOffset)
			tempSeeds = extendToTheRightDev(nodes[nodeIdx], read, readStart-(leftMatches-1), int(nodePos)-(leftMatches-1), false, tempSeeds)
			//log.Printf("After extendToTheRightDev rev:\n")
			//printSeedDev(tempSeeds)
			for _, tempSeed = range tempSeeds {
				//log.Printf("tempSeed.QueryStart = %d\n", tempSeed.QueryStart)
				finalSeeds = append(finalSeeds, extendToTheLeftDev(nodes[nodeIdx], read, tempSeed)...)
			}

		}
	}
	return finalSeeds
}
