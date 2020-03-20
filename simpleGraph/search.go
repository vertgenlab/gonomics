package simpleGraph

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
	//"fmt"
)

func AlignTraversalFwd(n *Node, seq []dna.Base, start int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, []uint32) {
	currentPath = append(currentPath, n.Id)
	var bestQueryEnd, queryEnd int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var path, bestPath []uint32

	if len(seq) >= ext {
		log.Fatalf("Error: the length of DNA sequence in previous nodes should not be enough to satisfy the desired extenion.\n")
	}
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	//log.Printf("len(seq)=%d, len(n.Seq)=%d, start=%d, targetLength=%d, basesToTake=%d\n", len(seq), len(n.Seq), start, targetLength, basesToTake)
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):targetLength], n.Seq[start:start+basesToTake])

	if availableBases >= ext || len(n.Next) == 0 {
		score, alignment, _, _, _, queryEnd = RightLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, queryEnd, currentPath
	} else {
		bestScore = -1
		tmpPath := make([]uint32, len(currentPath))
		copy(tmpPath, currentPath)
		for _, i := range n.Next {
			alignment, score, queryEnd, path = AlignTraversalFwd(i.Dest, s, 0, tmpPath, ext, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestQueryEnd = queryEnd
				bestPath = path
			}
		}
		return bestAlignment, bestScore, bestQueryEnd, bestPath
	}
}

/*
func AlignForwardTraversal(n *Node, seq []dna.Base, start int, seed *SeedDev, ext int, read *fastq.Fastq, m [][]int64, trace [][]rune) ([]*cigar.Cigar, *SeedDev) {

	var bestQueryEnd, queryEnd int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	//var path, bestPath []uint32

	if len(seq) >= ext {
		log.Fatalf("Error: the length of DNA sequence in previous nodes should not be enough to satisfy the desired extenion.\n")
	}
	var availableBases int = len(seq) + len(n.Seq) - start
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	//log.Printf("len(seq)=%d, len(n.Seq)=%d, start=%d, targetLength=%d, basesToTake=%d\n", len(seq), len(n.Seq), start, targetLength, basesToTake)
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):targetLength], n.Seq[start:start+basesToTake])

	if availableBases >= ext || len(n.Next) == 0 {

		score, alignment, _, _, _, queryEnd = RightLocal(s, read.Seq, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, seed
	} else {
		bestScore = -1
		newSeed := &SeedDev{TargetId: seed.TargetId, TargetStart: uint32(start), QueryStart: uint32(newQStart), Length: 0, PosStrand: seed.PosStrand, Next: nil, Prev: nil}
		for _, i := range n.Next {

			if availableBases < len(n.Seq) || availableBases < len(read.Seq[start:]) {

			}


			//tmpPath := make([]uint32, len(currentPath)) //TODO: should not need to copy path N times, but N-1 times
			copy(tmpPath, currentPath)
			alignment, score, queryEnd, path = AlignTraversalFwd(i.Dest, s, 0, tmpPath, ext, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestQueryEnd = queryEnd
				bestPath = path
			}
		}
		return bestAlignment, bestScore, bestQueryEnd, bestPath
	}
}*/

func AlignReverseGraphTraversal(n *Node, seq []dna.Base, refEnd int, currentPath []uint32, ext int, read []dna.Base, m [][]int64, trace [][]rune) ([]*cigar.Cigar, int64, int, int, []uint32) {
	currentPath = append([]uint32{n.Id}, currentPath...)
	var bestQueryStart, queryStart, refStart, bestRefStart int
	var bestScore, score int64
	var bestAlignment, alignment []*cigar.Cigar
	var path, bestPath []uint32

	var availableBases int = len(seq) + refEnd
	var targetLength int = common.Min(availableBases, ext)
	var basesToTake int = targetLength - len(seq)
	var s []dna.Base = make([]dna.Base, targetLength)
	copy(s[0:basesToTake], n.Seq[refEnd-basesToTake:refEnd])
	copy(s[basesToTake:targetLength], seq)

	//log.Printf("left(reverse) alignment: seq1=%s, seq2=%s\n", dna.BasesToString(s), dna.BasesToString(read))
	if availableBases >= ext || len(n.Next) == 0 {
		//log.Printf("at leaf, about to align, path is:%v\n", currentPath)
		score, alignment, refStart, _, queryStart, _ = LeftLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return alignment, score, refEnd - basesToTake + refStart, queryStart, currentPath
	} else {
		bestScore = -1
		tmp := make([]uint32, len(currentPath))
		copy(tmp, currentPath)
		for _, i := range n.Prev {

			alignment, score, refStart, queryStart, path = AlignReverseGraphTraversal(i.Dest, s, len(i.Dest.Seq), currentPath, ext, read, m, trace)
			if score > bestScore {
				bestScore = score
				bestAlignment = alignment
				bestRefStart = refStart
				bestQueryStart = queryStart
				bestPath = path
			}
		}
		return bestAlignment, bestScore, bestRefStart, bestQueryStart, bestPath
	}
}

func getSeqTraversal(curr *Node, seq []dna.Base, start int, extension int) [][]dna.Base {
	var answer [][]dna.Base
	if len(seq) >= extension {
		log.Fatalf("Error: the length of DNA sequence in previous nodes should not be enough to satisfy the desired extension.\n")
	}
	var availableBases int = len(curr.Seq) - start + len(seq)
	var targetLength int = common.Min(availableBases, extension)
	var basesToTake int = targetLength - len(seq)
	//log.Printf("len(seq)=%d, len(n.Seq)=%d, start=%d, targetLength=%d, basesToTake=%d\n", len(seq), len(curr.Seq), start, targetLength, basesToTake)
	var s []dna.Base = make([]dna.Base, targetLength)
	copy(s[0:len(seq)], seq)
	copy(s[len(seq):targetLength], curr.Seq[start:start+basesToTake])
	if availableBases >= extension || len(curr.Next) == 0 {
		if dna.CountBaseInterval(s, dna.N, 0, len(s)) == 0 {
			answer = append(answer, s)
		}

		//score, alignment, _, _, _, queryEnd = RightLocal(s, read, HumanChimpTwoScoreMatrix, -600, m, trace)
		return answer
	} else {

		for _, i := range curr.Next {
			answer = append(answer, getSeqTraversal(i.Dest, s, 0, extension)...)
		}
		return answer
	}
}

func lookingForSeeds(seedHash map[uint64][]*SeedBed, read *fastq.Fastq, seedLen int, stepSize int, posStrand bool, scoreMatrix [][]int64, gg *SimpleGraph) []*SeedDev {
	var codedSeq uint64 = 0
	var hits []*SeedDev = make([]*SeedDev, 0)
	var currSeed *SeedDev
	var currScore, bestScore int64 = 0, -1
	for subSeqStart := 0; subSeqStart < len(read.Seq)-seedLen+1; subSeqStart++ {
		if dna.CountBaseInterval(read.Seq, dna.N, subSeqStart, subSeqStart+seedLen) == 0 {
			codedSeq = dnaToNumber(read.Seq, subSeqStart, subSeqStart+seedLen)
			currHits := seedHash[codedSeq]
			for _, value := range currHits {
				currSeed = seedBedToSeedDev(value, uint32(subSeqStart), posStrand)
				extendSeedDev(currSeed, gg, read)
				currScore = BlastSeed(currSeed, read, scoreMatrix)
				if currScore > bestScore {
					hits = append(hits, currSeed)
				}
			}
		}
	}
	//log.Printf("Total of %d hits.\n", len(hits))
	return hits
}

/*
func devIndexGraph(genome *SimpleGraph, seedLen int, seedStep int) map[uint64][]*SeedBed {
	if seedLen < 2 || seedLen > 32 {
		log.Fatalf("Error: seed length needs to be greater than 1 and less than 33.  Got: %d\n", seedLen)
	}
	answer := make(map[uint64][]*SeedBed)
	var seqCode uint64
	var nodeIdx, pos int
	//var extension []*SeedBed
	for nodeIdx = 0; nodeIdx < len(genome.Nodes); nodeIdx++ {
		//log.Printf("Indexing node %d", nodeIdx)

		for pos = 0; pos < len(genome.Nodes[nodeIdx].Seq)-seedLen+1; pos += seedStep {
			if len(genome.Nodes[nodeIdx].Seq) < pos+seedLen {
				var extendSeq []dna.Base
				curr := SeedBed{Id: genome.Nodes[nodeIdx].Id, Start: uint32(pos), End: uint32(pos + len(genome.Nodes[nodeIdx].Seq)), Next: nil}
				//curr := SeedBed{Id: genome.Nodes[nodeIdx].Id, Start: uint32(pos), End: uint32(common.Min(len(genome.Nodes[nodeIdx].Seq), pos+seedLen)), Next: nil}
				_, answer = SeedBedGraph(genome, seedLen, &curr, extendSeq)

				//log.Printf("node after traversal %d", nodeIdx)

				//extension = SeedBedGraph(genome, seedLen, &curr, []dna.Base{})
				//fmt.Printf("extension seed len=%d", len(extension))
				//var extendSeq []dna.Base
				//for i = 0 ; i < len(extension);i++ {
				//	extendSeq = seedBedToSeq(extension[i], seedLen, genome)
				//	seqCode = dnaToNumber(extendSeq, 0, seedLen)
				//	answer[seqCode] = append(answer[seqCode], extension[i])
				//	nodeIdx = int(BedTail(extension[i]).Id)
				//}
			} else {
				if dna.CountBaseInterval(genome.Nodes[nodeIdx].Seq, dna.N, pos, pos+seedLen) == 0 {
					seqCode = dnaToNumber(genome.Nodes[nodeIdx].Seq, pos, pos+seedLen)
					curr := SeedBed{Id: genome.Nodes[nodeIdx].Id, Start: uint32(pos), End: uint32(pos + seedLen), Next: nil}
					answer[seqCode] = append(answer[seqCode], &curr)
				}
			}
		}
	}
	return answer
}

func BedTail(a *SeedBed) *SeedBed {
	if a.Next == nil {
		return a
	} else {
		return BedTail(a.Next)
	}
}

func seedBedToSeq(head *SeedBed, length int, gg *SimpleGraph) []dna.Base {
	var answer []dna.Base
	answer = append(answer, gg.Nodes[head.Id].Seq[:common.Min(length, len(gg.Nodes[head.Id].Seq))]...)
	if len(answer) < length && head.Next != nil {
		answer = append(answer, seedBedToSeq(head.Next, length, gg)...)
	}
	if len(answer) > length {
		answer = answer[:length]
	}
	//log.Printf("Checking to make sure len of sequence %d equals seedlen %d", len(answer), length)

	return answer
}

func SeedBedGraph(genome *SimpleGraph, seedLen int, graphSeed *SeedBed, leftOverSeq []dna.Base) (int, map[uint64][]*SeedBed) {
	if seedLen < 2 || seedLen > 32 {
		log.Fatalf("Error: seed length needs to be greater than 1 and less than 33.  Got: %d\n", seedLen)
	}
	var nodeIdx = int(graphSeed.Id)
	answer := make(map[uint64][]*SeedBed)
	var seqCode uint64
	var e int
	leftOverSeq = append(leftOverSeq, genome.Nodes[graphSeed.Id].Seq[graphSeed.Start:int(graphSeed.Start)+common.Min(seedLen-len(leftOverSeq), len(genome.Nodes[graphSeed.Id].Seq))]...)
	graphSeed.End = graphSeed.Start + uint32(len(leftOverSeq))
	if len(leftOverSeq) < seedLen {

		if len(genome.Nodes[graphSeed.Id].Next) > 0 {
			//var numBasesNeed int

			graphSeed.End = uint32(len(genome.Nodes[graphSeed.Id].Seq))
			//nextSeeds := make(map[uint64][]*SeedBed)
			for _, next := range genome.Nodes[graphSeed.Id].Next {
				nodeIdx = int(next.Dest.Id)
				nextGraphSeed := &SeedBed{Id: next.Dest.Id, Start: 0, End: uint32(common.Min(seedLen-len(leftOverSeq), len(next.Dest.Seq))), Next: nil}
				nodeIdx, answer = SeedBedGraph(genome, seedLen, nextGraphSeed, leftOverSeq)
				for edge := range answer {
					var tmp []*SeedBed
					for e = 0; e < len(answer[edge]); e++ {
						currSeed := &SeedBed{Id: graphSeed.Id, Start: graphSeed.Start, End: graphSeed.End, Next: answer[edge][e]}
						tmp = append(tmp, currSeed)
					}
					answer[edge] = tmp
				}
			}
		}
	} else {
		if dna.CountBaseInterval(leftOverSeq, dna.N, 0, len(leftOverSeq)) == 0 {
			//fmt.Printf("seqToCode: %s\n", dna.BasesToString(leftOverSeq))
			seqCode = dnaToNumber(leftOverSeq, 0, len(leftOverSeq))

			answer[seqCode] = append(answer[seqCode], graphSeed)

			//answer = append(answer, graphSeed)
		}
	}
	return nodeIdx, answer
}*/

func PointToBases(currSeq []dna.Base, incomingSeq []dna.Base, currStart int, incomingStart int, currEnd int, incomingEnd int, currFront bool) {
	newSize := currEnd + incomingEnd - currStart - incomingStart
	var i int
	if len(currSeq) < newSize {
		currSeq = append(currSeq, make([]dna.Base, newSize-len(currSeq))...)
	}
	if currFront {
		for i = currEnd; i < len(incomingSeq[incomingStart:incomingEnd]); i++ {
			currSeq[i] = incomingSeq[incomingStart+i]
		}
	}
	if !currFront {
		for i = 0; i < len(currSeq); i++ {
			currSeq[incomingEnd-incomingStart-1] = currSeq[i]
		}
		for i = 0; i < len(incomingSeq[incomingStart:incomingEnd]); i++ {
			currSeq[i] = incomingSeq[incomingStart+i]
		}
	}
}
