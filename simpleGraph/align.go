package simpleGraph

import("github.com/vertgenlab/gonomics/dna")

func indexGenome(genome []*Node, seedLen int) map[uint64][]uint64 {
	answer := make(map[uint64][]uint64)
	for chromIdx := 0; chromIdx < len(genome); chromIdx++ {
		for pos := 0; pos < len(genome[chromIdx].Seq)-seedLen+1; pos++ {
			seqCode := dnaToNumber(genome[chromIdx].Seq, pos, pos+seedLen)
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

func dnaToNumber(dna []dna.Base, start int, end int) uint64 {
	var answer uint64 = uint64(dna[start])
	for i := start + 1; i < end; i++ {
		answer = answer << 2
		answer = answer | uint64(dna[i])
	}
	return answer
}

func numberToChromAndPos(code uint64) (int, int) {
	var rightSideOnes uint64 = 4294967295
	var leftSideOnes uint64 = rightSideOnes << 32
	var chromIdx uint64 = code & leftSideOnes
	chromIdx = chromIdx >> 32
	var pos uint64 = code & rightSideOnes
	return int(chromIdx), int(pos)
}

