package genomeGraph

import (
	"fmt"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/numbers"
)

func createSeed(id uint32, readStart, nodeStart, length int, posStrand bool) *SeedDev {
	return &SeedDev{
		TargetId:    id,
		TargetStart: uint32(nodeStart),
		QueryStart:  uint32(readStart),
		Length:      uint16(length),
		PosStrand:   posStrand,
		TotalLength: uint16(length),
	}
}

// TODO: what about neg strand?
func perfectMatchBig(read fastq.FastqBig, scoreMatrix [][]int64) int64 {
	var fwdScore, revScore int64 = 0, 0
	for i, j := 0, len(read.Seq)-1; i <= j; i, j = i+1, j-1 {
		fwdScore = scoreMatrix[read.Seq[i]][read.Seq[i]]
		revScore = scoreMatrix[dna.ComplementSingleBase(read.Seq[j])][dna.ComplementSingleBase(read.Seq[j])]
	}
	return numbers.Max(fwdScore, revScore)
}

func scoreSeedSeq(seq []dna.Base, start uint32, end uint32, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for i := start; i < end; i++ {
		score += scoreMatrix[seq[i]][seq[i]]
	}
	return score
}

func scoreSeedPart(seed *SeedDev, read fastq.Fastq, scoreMatrix [][]int64) int64 {
	var score int64 = 0
	for _, base := range read.Seq[int(seed.QueryStart) : int(seed.QueryStart)+int(seed.Length)] {
		score += scoreMatrix[base][base]
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
