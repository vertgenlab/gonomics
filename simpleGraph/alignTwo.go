package simpleGraph

/*
import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
)

func GraphSmithWatermanMemPool(gg *SimpleGraph, read *fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, m [][]int64, trace [][]rune, memoryPool **SeedDev) *sam.SamAln {
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: "BZ:i:0"}
	//var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	//var i int
	var minTarget int
	//var minQuery int
	//var leftScore, rightScore int64
	var bestScore int64
	//var leftPath, rightPath []uint32

	var currScore int64 = 0
	perfectScore := perfectMatchBig(read)
	//extension := int(perfectScore/600) + len(read.Seq)

	var seeds *SeedDev
	seeds = findSeedsInSmallMapWithMemPool(seedHash, gg.Nodes, read, seedLen, perfectScore, memoryPool)
	SortSeedDevListByTotalLen(&seeds)

	var tailSeed *SeedDev
	var seedScore int64
	var currSeq []dna.Base
	var currSeed *SeedDev

	for currSeed = seeds; currSeed != nil && seedCouldBeBetter(int64(currSeed.TotalLength), bestScore, perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); currSeed = currSeed.Next {
		//tailSeed = toTail(currSeed)
		if currSeed.PosStrand {
			currSeq = read.Seq
		} else {
			currSeq = read.SeqRc
		}
		//leftAlignment, leftScore, minTarget, minQuery, leftPath = AlignReverseGraphTraversal(gg.Nodes[currSeed.TargetId], []dna.Base{}, int(currSeed.TargetStart), []uint32{}, extension, currSeq[:currSeed.QueryStart], m, trace)
		seedScore = scoreSeedSeq(currSeed, currSeq)
		//rightAlignment, rightScore, _, rightPath = AlignTraversalFwd(gg.Nodes[tailSeed.TargetId], []dna.Base{}, int(tailSeed.TargetStart+tailSeed.Length), []uint32{}, extension, currSeq[tailSeed.QueryStart+tailSeed.Length:], m, trace)

		//currScore = leftScore + seedScore + rightScore
		currScore = seedScore
		if currScore > bestScore {
			bestScore = currScore
			if currSeed.PosStrand {
				currBest.Flag = 0
			} else {
				currBest.Flag = 16
			}
			currBest.RName = fmt.Sprintf("%s_%d", gg.Nodes[currSeed.TargetId].Name, currSeed.TargetId)
			currBest.Pos = int64(minTarget) + 1
			//currBest.Cigar = cigar.CatCigar(cigar.AddCigar(leftAlignment, &cigar.Cigar{RunLength: int64(currSeed.TotalLength), Op: 'M'}), rightAlignment)
			//currBest.Cigar = AddSClip(minQuery, len(currSeq), currBest.Cigar)
			//currBest.Extra = "BZ:i:" + fmt.Sprint(bestScore) + "\tGP:Z:" + PathToString(CatPaths(CatPaths(leftPath, getSeedPath(currSeed)), rightPath), gg)
		}
	}
	if bestScore < 1200 {
		currBest.Flag = 4
	}
	tailSeed = toTail(seeds)
	tailSeed.Next = *memoryPool
	*memoryPool = seeds
	return &currBest
}

func GraphSmithWatermanBasic(gg *SimpleGraph, read *fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, m [][]int64, trace [][]rune) *sam.SamAln {
	var currBest sam.SamAln = sam.SamAln{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: "BZ:i:0"}
	var leftAlignment, rightAlignment []*cigar.Cigar = []*cigar.Cigar{}, []*cigar.Cigar{}
	var i, minTarget int
	var minQuery int
	var leftScore, rightScore, bestScore int64
	var leftPath, rightPath []uint32

	var currScore int64 = 0
	perfectScore := perfectMatchBig(read)
	extension := int(perfectScore/600) + len(read.Seq)

	var seeds []*SeedDev
	seeds = findSeedsInSmallMap(seedHash, gg.Nodes, read, seedLen, perfectScore)
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
}

//TODO: what about neg strand?
func perfectMatchBig(read *fastq.FastqBig) int64 {
	var perfectScore int64 = 0
	for i := 0; i < len(read.Seq); i++ {
		perfectScore += HumanChimpTwoScoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return perfectScore
}

func scoreSeedSeq(seed *SeedDev, seq []dna.Base) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		score += HumanChimpTwoScoreMatrix[seq[i]][seq[i]]
	}
	return score
}

func scoreSeedFastqBig(seed *SeedDev, read *fastq.FastqBig) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		if seed.PosStrand {
			score += HumanChimpTwoScoreMatrix[read.Seq[i]][read.Seq[i]]
		} else {
			score += HumanChimpTwoScoreMatrix[read.SeqRc[i]][read.SeqRc[i]]
		}
	}
	return score
}

func scoreSeedPart(seed *SeedDev, read *fastq.Fastq) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		score += HumanChimpTwoScoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return score
}

func scoreSeed(seed *SeedDev, read *fastq.Fastq) int64 {
	var score int64 = 0
	for ; seed != nil; seed = seed.NextPart {
		score += scoreSeedPart(seed, read)
	}
	return score
}*/
