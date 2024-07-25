package genomeGraph

import (
	//	"github.com/vertgenlab/gonomics/fastq"
	"log"

	"github.com/vertgenlab/gonomics/dna"
)

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

func IndexGenomeIntoMap(genome []Node, seedLen int, seedStep int) map[uint64][]uint64 {
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

func indexGenomeIntoSlice(genome []Node, seedLen int, seedStep int) [][]uint64 {
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
