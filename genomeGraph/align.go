package genomeGraph

import (
	"fmt"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
)

func GraphSmithWatermanMemPool(gg *GenomeGraph, read fastq.FastqBig, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, m [][]int64, trace [][]rune, memoryPool **SeedDev) sam.Sam {
	var currBest sam.Sam = sam.Sam{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []cigar.Cigar{{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: "", Extra: "BZ:i:0\tGP:Z:-1"}
	var leftAlignment, rightAlignment []cigar.Cigar
	var minTarget int
	var minQuery int
	var leftScore, rightScore int64
	var bestScore int64
	var leftPath, rightPath, bestPath []uint32
	var currScore int64 = 0
	perfectScore := perfectMatchBig(read, scoreMatrix)

	var seeds []*SeedDev
	seeds = findSeedsInSmallMapWithMemPool(seedHash, gg.Nodes, read, seedLen, perfectScore, scoreMatrix)
	SortSeedDevByLen(seeds)
	var tailSeed *SeedDev
	var seedScore int64
	var currSeq []dna.Base
	var currSeed *SeedDev

	for i := 0; i < len(seeds) && seedCouldBeBetter(int64(seeds[i].TotalLength), bestScore, perfectScore, int64(len(read.Seq)), 100, 90, -196, -296); i++ {
		currSeed = seeds[i]
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
			currBest.RName = fmt.Sprintf("%d", bestPath[0])
			currBest.Pos = uint32(minTarget + 1)
			currBest.Extra = "BZ:i:" + fmt.Sprint(bestScore) + "\tGP:Z:" + PathToString(CatPaths(CatPaths(leftPath, getSeedPath(currSeed)), rightPath))

			currBest.Cigar = cigar.CatCigar(cigar.AddCigar(leftAlignment, cigar.Cigar{RunLength: int(currSeed.TotalLength), Op: 'M'}), rightAlignment)
			currBest.Cigar = AddSClip(minQuery, len(currSeq), currBest.Cigar)
		}
	}
	if bestScore < 1200 {
		currBest.Flag = 4
	}
	return currBest
}

// TODO: what about neg strand?
func perfectMatchBig(read fastq.FastqBig, scoreMatrix [][]int64) int64 {
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

func scoreSeedFastqBig(seed *SeedDev, read fastq.FastqBig, scoreMatrix [][]int64) int64 {
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

func scoreSeedPart(seed *SeedDev, read fastq.Fastq, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for i := seed.QueryStart; i < seed.QueryStart+seed.Length; i++ {
		score += scoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return score
}

func scoreSeed(seed *SeedDev, read fastq.Fastq, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for ; seed != nil; seed = seed.NextPart {
		score += scoreSeedPart(seed, read, scoreMatrix)
	}
	return score
}

var HumanChimpTwoScoreMatrix = [][]int64{
	{90, -330, -236, -356, -208},
	{-330, 100, -318, -236, -196},
	{-236, -318, 100, -330, -196},
	{-356, -236, -330, 90, -208},
	{-208, -196, -196, -208, -202},
}

func AddSClip(front int, lengthOfRead int, cig []cigar.Cigar) []cigar.Cigar {
	var runLen int = cigar.QueryLength(cig)
	if runLen < lengthOfRead {
		answer := make([]cigar.Cigar, 0, len(cig)+2)
		if front > 0 {
			answer = append(answer, cigar.Cigar{RunLength: front, Op: 'S'})
		}
		answer = append(answer, cig...)
		if front+int(cigar.QueryLength(cig)) < lengthOfRead {
			answer = append(answer, cigar.Cigar{RunLength: lengthOfRead - front - runLen, Op: 'S'})
		}
		return answer
	} else {
		return cig
	}
}

// perfect match
func perfectMatch(read fastq.Fastq, scoreMatrix [][]int64) int64 {
	var perfectScore int64 = 0
	for i := 0; i < len(read.Seq); i++ {
		perfectScore += scoreMatrix[read.Seq[i]][read.Seq[i]]
	}
	return perfectScore
}

/*func NodesHeader(ref []*Node) *sam.Header {
	var header sam.Header
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string
	for i := 0; i < len(ref); i++ {
		words = fmt.Sprintf("@SQ\tSN:%s_%d\tLN:%d", ref[i].Name, ref[i].Id, len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, &chromInfo.ChromInfo{Name: ref[i].Name, Size: len(ref[i].Seq)})
	}
	return &header
}*/

func ChromAndPosToNumber(chrom int, start int) uint64 {
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

func ViewMatrix(m [][]int64) string {
	var message string = ""
	message += fmt.Sprintf("\t\t %d\t%d\t%d\t%d\n\t\t%d\t%d\t%d\t%d\n\t\t%d\t%d\t%d\t%d\n\t\t%d\t%d\t%d\t %d\n", m[0][0], m[0][1], m[0][2], m[0][3], m[1][0], m[1][1], m[1][2], m[1][3], m[2][0], m[2][1], m[2][2], m[2][3], m[3][0], m[3][1], m[3][2], m[3][3])
	return message
}
