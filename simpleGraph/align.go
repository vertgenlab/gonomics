package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
)

//Uses small mem pool
func GraphSmithWaterman(gg *SimpleGraph, read *fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]rune, memoryPool **SeedDev) *sam.SamAln {
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: "BZ:i:0\tGP:Z:-1"}
	var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	var minTarget int
	var minQuery int
	var leftScore, rightScore int64 = 0, 0
	var bestScore int64
	var leftPath, rightPath, bestPath []uint32
	var currScore int64 = 0
	perfectScore := perfectMatchBig(read, scoreMatrix)
	extension := int(perfectScore/600) + len(read.Seq)
	var seeds *SeedDev
	seeds = findSeedsInSmallMapWithMemPool(seedHash, gg.Nodes, read, seedLen, perfectScore, scoreMatrix, memoryPool)
	SortSeedDevListByTotalLen(&seeds)
	var tailSeed *SeedDev
	var seedScore int64
	var currSeq []dna.Base
	var currSeed *SeedDev
	for currSeed = seeds; currSeed != nil && seedCouldBeBetterScores(int64(currSeed.TotalLength), bestScore, perfectScore, int64(len(read.Seq)), scoreMatrix); currSeed = currSeed.Next {
		tailSeed = getLastPart(currSeed)
		if currSeed.PosStrand {
			currSeq = read.Seq
		} else {
			currSeq = read.SeqRc
		}
		if int(currSeed.TotalLength) == len(currSeq) {
			currScore = seedScore
			minTarget = int(currSeed.TargetStart)
			minQuery = int(currSeed.QueryStart)
		} else {
			leftAlignment, leftScore, minTarget, minQuery, leftPath = AlignReverseGraphTraversal(gg.Nodes[currSeed.TargetId], []dna.Base{}, int(currSeed.TargetStart), []uint32{}, extension-int(currSeed.TotalLength), currSeq[:currSeed.QueryStart], m, trace)
			rightAlignment, rightScore, _, _, rightPath = AlignTraversalFwd(gg.Nodes[tailSeed.TargetId], []dna.Base{}, int(tailSeed.TargetStart+tailSeed.Length), []uint32{}, extension-int(currSeed.TotalLength), currSeq[tailSeed.QueryStart+tailSeed.Length:], m, trace)
		}
		seedScore = scoreSeedSeq(currSeq, currSeed.QueryStart, tailSeed.QueryStart+tailSeed.Length, scoreMatrix)
		currScore = leftScore + seedScore + rightScore
		if currScore > bestScore {
			bestPath = CatPaths(CatPaths(leftPath, getSeedPath(currSeed)), rightPath)
			bestScore = currScore
			if currSeed.PosStrand {
				currBest.Flag = 0
			} else {
				currBest.Flag = 16
			}
			currBest.Seq = currSeq
			currBest.Qual = string(read.Qual)
			currBest.RName = gg.Nodes[bestPath[0]].Name
			currBest.Pos = int64(minTarget) + 1
			if gg.Nodes[bestPath[0]].Info != nil {
				currBest.Extra = fmt.Sprintf("BZ:i:%d\tGP:Z:%s\tXO:i:%d", bestScore, PathToString(bestPath), gg.Nodes[bestPath[0]].Info.Start-1)
			} else {
				currBest.Extra = fmt.Sprintf("BZ:i:%d\tGP:Z:%s", bestScore, PathToString(bestPath))
			}
			currBest.Cigar = AddSClip(minQuery, len(currSeq), cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(sumLen(currSeed)), Op: 'M'}), rightAlignment))
		}
	}
	if bestScore < 1200 {
		currBest.Flag = 4
	}
	if seeds != nil {
		tailSeed = toTail(seeds)
		tailSeed.Next = *memoryPool
		*memoryPool = seeds
	}
	return &currBest
}

//let me know if i figured this out right...
//returns maxMatch, minMatch, leastSevereMismatch, leastSevereMatchMismatchChange
//100, 90, -196, -296
func calcSeedMismatch(scoreMatrix [][]int64) (int64, int64, int64, int64) {
	//returns maxMatch, minMatch, leastSevereMismatch, leastSevereMatchMismatchChange
	//100, 90, -196, -296
	//(2,2), (1,1), (1,4), (1,4)-(2,2)
	return scoreMatrix[2][2], scoreMatrix[1][1], scoreMatrix[1][4], (scoreMatrix[1][4] - scoreMatrix[2][2])
}

func seedCouldBeBetterScores(seedLen int64, currBestScore int64, perfectScore int64, queryLen int64, scoreMatrix [][]int64) bool {
	maxMatch, minMatch, leastSevereMismatch, leastSevereMatchMismatchChange := calcSeedMismatch(scoreMatrix)
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

func GraphSmithWatermanMemPool(gg *SimpleGraph, read *fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]rune, memoryPool **SeedDev) *sam.SamAln {
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: "", Extra: "BZ:i:0\tGP:Z:-1"}
	var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	var minTarget int
	var minQuery int
	var leftScore, rightScore int64
	var bestScore int64
	var leftPath, rightPath, bestPath []uint32
	var currScore int64 = 0
	perfectScore := perfectMatchBig(read, scoreMatrix)
	extension := int(perfectScore/600) + len(read.Seq)
	var seeds *SeedDev
	seeds = findSeedsInSmallMapWithMemPool(seedHash, gg.Nodes, read, seedLen, perfectScore, scoreMatrix, memoryPool)
	SortSeedDevListByTotalLen(&seeds)
	var tailSeed *SeedDev
	var seedScore int64
	var currSeq []dna.Base
	var currSeed *SeedDev
	for currSeed = seeds; currSeed != nil && seedCouldBeBetter(int64(currSeed.TotalLength), bestScore, perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); currSeed = currSeed.Next {
		tailSeed = getLastPart(currSeed)
		if currSeed.PosStrand {
			currSeq = read.Seq
		} else {
			currSeq = read.SeqRc
		}
		if int(currSeed.TotalLength) == len(currSeq) {
			currScore = seedScore
			leftScore = 0
			minTarget = int(currSeed.TargetStart)
			minQuery = int(currSeed.QueryStart)
			rightScore = 0
		} else {
			leftAlignment, leftScore, minTarget, minQuery, leftPath = AlignReverseGraphTraversal(gg.Nodes[currSeed.TargetId], []dna.Base{}, int(currSeed.TargetStart), []uint32{}, extension-int(currSeed.TotalLength), currSeq[:currSeed.QueryStart], m, trace)
			rightAlignment, rightScore, _, _, rightPath = AlignTraversalFwd(gg.Nodes[tailSeed.TargetId], []dna.Base{}, int(tailSeed.TargetStart+tailSeed.Length), []uint32{}, extension-int(tailSeed.TotalLength), currSeq[tailSeed.QueryStart+tailSeed.Length:], m, trace)
		}
		seedScore = scoreSeedSeq(currSeq, currSeed.QueryStart, tailSeed.QueryStart+tailSeed.Length, scoreMatrix)
		currScore = leftScore + seedScore + rightScore
		if currScore > bestScore {
			bestPath = CatPaths(CatPaths(leftPath, getSeedPath(currSeed)), rightPath)
			bestScore = currScore
			if currSeed.PosStrand {
				currBest.Flag = 0
			} else {
				currBest.Flag = 16
			}
			currBest.Seq = currSeq // unsure why this line was lost
			currBest.Qual = string(read.Qual)
			currBest.RName = gg.Nodes[bestPath[0]].Name
			currBest.Pos = int64(minTarget) + 1
			currBest.Extra = "BZ:i:" + fmt.Sprint(bestScore) + "\tGP:Z:" + PathToString(CatPaths(CatPaths(leftPath, getSeedPath(currSeed)), rightPath))
			if gg.Nodes[bestPath[0]].Info != nil {
				currBest.Extra += fmt.Sprintf("\tXO:i:%d", gg.Nodes[bestPath[0]].Info.Start-1)
				//currBest.Pos += int64(gg.Nodes[bestPath[0]].Info.Start) - 1
			}
			currBest.Cigar = cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(sumLen(currSeed)), Op: 'M'}), rightAlignment)
			currBest.Cigar = AddSClip(minQuery, len(currSeq), currBest.Cigar)
		}
	}
	if bestScore < 1200 {
		currBest.Flag = 4
	}
	if seeds != nil {
		tailSeed = toTail(seeds)
		tailSeed.Next = *memoryPool
		*memoryPool = seeds
	}
	return &currBest
}

/*func GraphSmithWatermanBasic(gg *SimpleGraph, read *fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]rune) *sam.SamAln {
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: make([]dna.Base, 0, len(read.Seq)), Qual: "", Extra: "BZ:i:0\tGP:Z:-1"}
	var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	var i, minTarget int
	var minQuery int
	var leftScore, rightScore, bestScore int64
	var leftPath, rightPath []uint32
	var currScore int64 = 0
	perfectScore := perfectMatchBig(read)
	extension := int(perfectScore/600) + len(read.Seq)
	var seeds []*SeedDev
	seeds = findSeedsInSmallMap(seedHash, gg.Nodes, read, seedLen, perfectScore, scoreMatrix)
	SortSeedDevByTotalLen(seeds)
	var tailSeed *SeedDev
	var seedScore int64
	var currSeq []dna.Base
	for i = 0; i < len(seeds) && seedCouldBeBetter(int64(seeds[i].TotalLength), bestScore, perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		tailSeed = toTail(seeds[i])
		if seeds[i].PosStrand {
			currSeq = read.Seq
		} else {
			currSeq = read.SeqRc
		}
		leftAlignment, leftScore, minTarget, minQuery, leftPath = AlignReverseGraphTraversal(gg.Nodes[seeds[i].TargetId], []dna.Base{}, int(seeds[i].TargetStart), []uint32{}, extension, currSeq[:seeds[i].QueryStart], m, trace)
		seedScore = scoreSeedSeq(seeds[i], currSeq)
		rightAlignment, rightScore, _, rightPath = AlignTraversalFwd(gg.Nodes[tailSeed.TargetId], []dna.Base{}, int(tailSeed.TargetStart+tailSeed.Length), []uint32{}, extension, currSeq[tailSeed.QueryStart+tailSeed.Length:], m, trace)
		currScore = leftScore + seedScore + rightScore
		if currScore > bestScore {
			bestScore = currScore
			if seeds[i].PosStrand {
				currBest.Flag = 0
			} else {
				currBest.Flag = 16
			}
			currBest.RName = fmt.Sprintf("%s_%d", gg.Nodes[seeds[i].TargetId].Name, seeds[i].TargetId)
			currBest.Pos = int64(minTarget) + 1
			currBest.Cigar = cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(seeds[i].TotalLength), Op: 'M'}), rightAlignment)
			currBest.Cigar = AddSClip(minQuery, len(currSeq), currBest.Cigar)
			currBest.Extra = "BZ:i:" + fmt.Sprint(bestScore) + "\tGP:Z:" + PathToString(CatPaths(CatPaths(leftPath, getSeedPath(seeds[i])), rightPath), gg)
		}
	}
	if bestScore < 1200 {
		currBest.Flag = 4
	}
	return &currBest
}*/

//TODO: what about neg strand?
func perfectMatchBig(read *fastq.FastqBig, scoreMatrix [][]int64) int64 {
	var perfectScore int64 = 0
	for i := 0; i < len(read.Seq); i++ {
		perfectScore += scoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return perfectScore
}

func scoreSeedSeq(seq []dna.Base, start uint32, end uint32, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for i := start; i < end; i++ {
		score += scoreMatrix[seq[i]][seq[i]]
	}
	return score
}

func scoreSeedFastqBig(seed *SeedDev, read *fastq.FastqBig, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		if seed.PosStrand {
			score += scoreMatrix[read.Seq[i]][read.Seq[i]]
		} else {
			score += scoreMatrix[read.SeqRc[i]][read.SeqRc[i]]
		}
	}
	return score
}

func scoreSeedPart(seed *SeedDev, read *fastq.Fastq, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		score += scoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return score
}

func scoreSeed(seed *SeedDev, read *fastq.Fastq, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for ; seed != nil; seed = seed.NextPart {
		score += scoreSeedPart(seed, read, scoreMatrix)
	}
	return score
}

func oldGraphSmithWaterman(gg *SimpleGraph, read *fastq.Fastq, seedHash map[uint64][]*SeedBed, seedLen int, stepSize int, m [][]int64, trace [][]rune) *sam.SamAln {
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: make([]dna.Base, 0, len(read.Seq)), Qual: "", Extra: "BZ:i:0\tGP:Z:-1"}
	var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	var i, minTarget int
	var minQuery int
	var leftScore, rightScore, bestScore int64
	var leftPath, rightPath, bestPath []uint32

	var currScore int64 = 0
	perfectScore := perfectMatch(read, HumanChimpTwoScoreMatrix)
	extension := int(perfectScore/600) + len(read.Seq)

	var currRead *fastq.Fastq = nil
	var seeds []*SeedDev = findSeedsInMapDev(seedHash, read, seedLen, stepSize, true)
	//var seeds []*SeedDev = lookingForSeeds(seedHash, read, seedLen, stepSize, true, HumanChimpTwoScoreMatrix, gg)
	seeds = GraphDictionary(seeds, gg, read)

	revCompRead := fastq.Copy(read)
	fastq.ReverseComplement(revCompRead)
	var revCompSeeds []*SeedDev = findSeedsInMapDev(seedHash, revCompRead, seedLen, stepSize, false)
	//var revCompSeeds []*SeedDev = lookingForSeeds(seedHash, revCompRead, seedLen, stepSize, true, HumanChimpTwoScoreMatrix, gg)
	revCompSeeds = GraphDictionary(revCompSeeds, gg, revCompRead)

	seeds = append(seeds, revCompSeeds...)
	SortSeedExtended(seeds)
	var tailSeed *SeedDev
	var seedScore int64

	for i = 0; i < len(seeds) && seedCouldBeBetter(int64(seeds[i].TotalLength), bestScore, perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		tailSeed = toTail(seeds[i])
		if seeds[i].PosStrand {
			currRead = read
		} else {
			currRead = revCompRead
		}

		leftAlignment, leftScore, minTarget, minQuery, leftPath = AlignReverseGraphTraversal(gg.Nodes[seeds[i].TargetId], []dna.Base{}, int(seeds[i].TargetStart), []uint32{}, extension, currRead.Seq[:seeds[i].QueryStart], m, trace)
		//log.Printf("NodeLen=%d, TargetStart=%d, length=%d\n", len(gg.Nodes[tailSeed.TargetId].Seq), tailSeed.TargetStart, tailSeed.Length)
		seedScore = BlastSeed(seeds[i], currRead, HumanChimpTwoScoreMatrix)
		rightAlignment, rightScore, _, _, rightPath = AlignTraversalFwd(gg.Nodes[tailSeed.TargetId], []dna.Base{}, int(tailSeed.TargetStart+tailSeed.Length), []uint32{}, extension, currRead.Seq[tailSeed.QueryStart+tailSeed.Length:], m, trace)
		currScore = leftScore + seedScore + rightScore

		if currScore > bestScore {
			bestPath = CatPaths(CatPaths(leftPath, getSeedPath(seeds[i])), rightPath)
			bestScore = currScore
			if seeds[i].PosStrand {
				currBest.Flag = 0
			} else {
				currBest.Flag = 16
			}
			currBest.Seq = currRead.Seq
			currBest.Qual = string(currRead.Qual)
			//if gg.Nodes[bestPath[0]].Info != nil {
			//	currBest.RName = fmt.Sprintf("%s.%d.%d", gg.Nodes[bestPath[0]].Name, gg.Nodes[bestPath[0]].Id, gg.Nodes[bestPath[0]].Info.Start)
			//}
			currBest.RName = gg.Nodes[bestPath[0]].Name
			//currBest.RName = fmt.Sprintf("%s_%d", gg.Nodes[bestPath[0]].Name, gg.Nodes[bestPath[0]].Id)
			currBest.Pos = int64(minTarget) + 1
			currBest.Extra = "BZ:i:" + fmt.Sprint(bestScore) + "\tGP:Z:" + PathToString(CatPaths(CatPaths(leftPath, getSeedPath(seeds[i])), rightPath))
			if gg.Nodes[bestPath[0]].Info != nil {
				currBest.Extra += fmt.Sprintf("\tXO:i:%d", gg.Nodes[bestPath[0]].Info.Start-1)
				currBest.Pos += int64(gg.Nodes[bestPath[0]].Info.Start) - 1
			}
			currBest.Cigar = cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(sumLen(seeds[i])), Op: 'M'}), rightAlignment)
			currBest.Cigar = AddSClip(minQuery, len(currRead.Seq), currBest.Cigar)

		}
	}
	if bestScore < 1200 {
		currBest.Flag = 4
	}

	return &currBest
}

var HumanChimpTwoScoreMatrix = [][]int64{
	{90, -330, -236, -356, -208},
	{-330, 100, -318, -236, -196},
	{-236, -318, 100, -330, -196},
	{-356, -236, -330, 90, -208},
	{-208, -196, -196, -208, -202},
}

func AddSClip(front int, lengthOfRead int, cig []*cigar.Cigar) []*cigar.Cigar {
	var runLen int64 = cigar.QueryLength(cig)
	if runLen < int64(lengthOfRead) {
		answer := make([]*cigar.Cigar, 0, len(cig)+2)
		if front > 0 {
			answer = append(answer, &cigar.Cigar{RunLength: int64(front), Op: 'S'})
		}
		answer = append(answer, cig...)
		if front+int(cigar.QueryLength(cig)) < lengthOfRead {
			answer = append(answer, &cigar.Cigar{RunLength: int64(lengthOfRead-front) - runLen, Op: 'S'})
		}
		return answer
	} else {
		return cig
	}
}

//perfect match
func perfectMatch(read *fastq.Fastq, scoreMatrix [][]int64) int64 {
	var perfectScore int64 = 0
	for i := 0; i < len(read.Seq); i++ {
		perfectScore += scoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return perfectScore
}

func NodesHeader(ref []*Node) *sam.SamHeader {
	var header sam.SamHeader
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string

	for i := 0; i < len(ref); i++ {
		//words = "@SQ\tSN:" + strconv.FormatInt(ref[i].Id, 10) + "\tLN:" + strconv.Itoa(len(ref[i].Seq))
		//words = "@SQ\tSN:" + ref[i].Name + "_" + fmt.Sprint(ref[i].Id) + "\tLN:" + fmt.Sprint(len(ref[i].Seq))
		words = fmt.Sprintf("@SQ\tSN:%s_%d\tLN:%d", ref[i].Name, ref[i].Id, len(ref[i].Seq))
		//words = fmt.Sprintf("@SQ\tSN:%s\tLN:%d", ref[i].Name, len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, &chromInfo.ChromInfo{Name: ref[i].Name, Size: int64(len(ref[i].Seq))})
	}
	return &header
}
func indexGenome(genome []*Node, seedLen int) map[uint64][]uint64 {
	answer := make(map[uint64][]uint64)
	var seqCode uint64
	var chromIdx, pos int
	for chromIdx = 0; chromIdx < len(genome); chromIdx++ {

		for pos = 0; pos < len(genome[chromIdx].Seq)-seedLen+1; pos++ {
			seqCode = dnaToNumber(genome[chromIdx].Seq, pos, pos+seedLen)
			answer[seqCode] = append(answer[seqCode], chromAndPosToNumber(chromIdx, pos))
		}
	}
	return answer
}

func chromAndPosToNumber(chrom int, start int) uint64 {
	var chromCode uint64 = uint64(chrom)
	chromCode = chromCode << 32
	var answer uint64 = chromCode | uint64(start)
	return answer
}

func dnaToNumber(seq []dna.Base, start int, end int) uint64 {
	var answer uint64 = uint64(seq[start])
	for i := start + 1; i < end; i++ {
		answer = answer << 2
		answer = answer | uint64(seq[i])
	}
	return answer
}

func numberToChromAndPos(code uint64) (int64, int64) {
	var rightSideOnes uint64 = 4294967295
	var leftSideOnes uint64 = rightSideOnes << 32
	var chromIdx uint64 = code & leftSideOnes
	chromIdx = chromIdx >> 32
	var pos uint64 = code & rightSideOnes
	return int64(chromIdx), int64(pos)
}

/*

func maskSeedLen(seq []dna.Base, start int, end int, numMismatch int) []uint64 {
	var answer []uint64
	var maskOne uint64 = 1
	var maskTwo uint64 = 2
	var maskThree uint64 = 3
	for i := 1; i < numMismatch; i++ {
		maskOne = maskOne << 1
		maskTwo = maskTwo << 2
		maskThree = maskThree << 3
	}
	answer = append(answer, dnaToNumber(seq, start, end)^maskOne)
	answer = append(answer, dnaToNumber(seq, start, end)^maskTwo)
	answer = append(answer, dnaToNumber(seq, start, end)^maskThree)
	return answer
}*/

/*
func MapSingleFastq(ref []*Node, chromPosHash map[uint64][]uint64, read *fastq.Fastq, seedLen int, m [][]int64, trace [][]rune) *sam.SamAln {
	//var seedBeds []*bed.Bed
	var seeds []Seed = make([]Seed, 256)
	var codedPositions []uint64
	var bStart, bEnd int64
	var maxScore int64
	var extension int64
	var chrom, pos int64
	var bases int
	var i int
	var subRead, hits int
	var score, bestScore int64 = 0, 0
	var lowRef, lowQuery, highQuery int64
	var alignment []*cigar.Cigar
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 0, RName: "", Pos: 0, MapQ: 255, RNext: "*", PNext: 0, TLen: 0, Seq: []dna.Base{}, Qual: "", Extra: ""}
	for bases = 0; bases < len(read.Seq); bases++ {
		maxScore += HumanChimpTwoScoreMatrix[read.Seq[bases]][read.Seq[bases]]
	}
	extension = int64(maxScore/600) + int64(len(read.Seq))
	for subRead = 0; subRead < len(read.Seq)-seedLen+1; subRead++ {
		codedPositions = chromPosHash[dnaToNumber(read.Seq, subRead, subRead+seedLen)]
		for hits = 0; hits < len(codedPositions); hits++ {
			chrom, pos = numberToChromAndPos(codedPositions[hits])
			bStart = pos - extension
			if bStart < 0 {
				bStart = 0
			}
			bEnd = pos + extension
			if bEnd > int64(len(ref[chrom].Seq)) {
				bEnd = int64(len(ref[chrom].Seq))
			}
			//seedBeds = append(seedBeds, &bed.Bed{Chrom: strconv.FormatInt(chrom, 10), ChromStart: bStart, ChromEnd: bEnd, Score: 1})
			seeds = addSeed(seeds, chrom, bStart, bEnd)
		}
	}
	//seedBeds = bed.MergeBeds(seedBeds)
	//log.Printf("Length of seed beds, %d", len(seedBeds))
	for i = 0; i < len(seeds); i++ {
		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[seeds[i].Id].Seq[seeds[i].Start:seeds[i].End], read.Seq, HumanChimpTwoScoreMatrix, -600, m, trace)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = fmt.Sprintf("%d", seeds[i].Id)
			currBest.Pos = lowRef + seeds[i].Start
			currBest.Cigar = SClipCigar(lowQuery, highQuery, int64(len(read.Seq)), alignment)
			currBest.Seq = read.Seq
			currBest.Qual = string(read.Qual)
		}
	}
	seeds = seeds[0:0]
	//seedBeds = [0:0]
	fastq.ReverseComplement(read)
	//extension = int64(maxScore/600) + int64(len(reverse.Seq))
	for subRead = 0; subRead < len(read.Seq)-seedLen+1; subRead++ {
		codedPositions = chromPosHash[dnaToNumber(read.Seq, subRead, subRead+seedLen)]
		for hits = 0; hits < len(codedPositions); hits++ {
			chrom, pos = numberToChromAndPos(codedPositions[hits])
			if chrom >= int64(len(ref)) {
				log.Printf("hashing chrom had index out of bounds b/c of Ns")
			} else {
				bStart = pos - extension
				if bStart < 0 {
					bStart = 0
				}
				bEnd = pos + extension
				if bEnd > int64(len(ref[chrom].Seq)) {
					bEnd = int64(len(ref[chrom].Seq))
				}
				//seedBeds = append(seedBeds, &bed.Bed{Chrom: strconv.FormatInt(chrom, 10), ChromStart: bStart, ChromEnd: bEnd, Score: 1})
				seeds = addSeed(seeds, chrom, bStart, bEnd)
			}
		}
	}
	//seedBeds = bed.MergeBeds(seedBeds)
	//log.Printf("Length of reverse seed beds, %d", len(seedBeds))
	for i = 0; i < len(seeds); i++ {
		score, alignment, lowRef, _, lowQuery, highQuery = SmithWaterman(ref[seeds[i].Id].Seq[seeds[i].Start:seeds[i].End], read.Seq, HumanChimpTwoScoreMatrix, -600, m, trace)
		if score > bestScore {
			bestScore = score
			currBest.Flag = 0
			currBest.RName = fmt.Sprintf("%d", seeds[i].Id)
			currBest.Pos = lowRef + seeds[i].Start
			currBest.Cigar = SClipCigar(lowQuery, highQuery, int64(len(read.Seq)), alignment)
			currBest.Seq = read.Seq
			currBest.Qual = string(read.Qual)
			//fmt.Printf("alignment:\n%s\n", cigar.LocalView(ref[seeds[i].Id].Seq[seeds[i].Start:seeds[i].End], read.Seq, alignment, highRef))
		}
	}
	//seeds = seeds[0:0]
	//log.Println(currBest.RName, "\t", currBest.Pos, "\t", read.Name, "\t", bestScore, "\t", cigar.ToString(currBest.Cigar))
	return &currBest
}*/
