package genomeGraph

import (
	//	"github.com/vertgenlab/gonomics/fastq"
	"log"
	"sort"
	"sync"

	"github.com/vertgenlab/gonomics/dna"
)

const (
	rightMask uint64 = 0xFFFFFFFF      // 32 ones, to mask the right side (position)
	leftMast  uint64 = rightMask << 32 // 32 ones shifted, to mask the left side (chromosome)
)

func IndexGenomeIntoMap(genome []Node, seedLen int, seedStep int) map[uint64][]uint64 {
	if seedLen < 2 || seedLen > 32 {
		log.Fatalf("Error: seed length needs to be greater than 1 and less than 33.  Got: %d\n", seedLen)
	}
	var wg sync.WaitGroup
	answer := make(map[uint64][]uint64)
	var mu sync.Mutex

	for nodeIdx := range genome {
		wg.Add(1)
		go func(nodeIdx int) {
			defer wg.Done()
			localMap := make(map[uint64][]uint64)
			var seqCode, locationCode uint64
			var pos int
			for pos = 0; pos < len(genome[nodeIdx].Seq)-seedLen+1; pos += seedStep {
				if dna.CountBaseInterval(genome[nodeIdx].Seq, dna.N, pos, pos+seedLen) == 0 {
					seqCode = dnaToNumber(genome[nodeIdx].Seq, pos, pos+seedLen)
					localMap[seqCode] = append(localMap[seqCode], ChromAndPosToNumber(nodeIdx, pos))
				}
			}
			for ; pos < len(genome[nodeIdx].Seq); pos += seedStep {
				locationCode = ChromAndPosToNumber(nodeIdx, pos)
				for edgeIdx := 0; edgeIdx < len(genome[nodeIdx].Next); edgeIdx++ {
					indexGenomeIntoMapHelper(genome[nodeIdx].Seq[pos:], genome[nodeIdx].Next[edgeIdx].Dest, locationCode, seedLen, localMap)
				}
			}

			// Synchronize map update
			mu.Lock()
			for k, v := range localMap {
				answer[k] = append(answer[k], v...)
			}
			mu.Unlock()
		}(nodeIdx)
	}

	wg.Wait()
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

// TODO: this does not take into account breaking up seeds by gaps instead of mismatches
// similar calculations could also be used as the parameters to a banded alignment.
func seedCouldBeBetter(seedLen int64, currBestScore int64, perfectScore int64, queryLen int64, config *GraphSettings) bool {
	seeds := queryLen / (seedLen + 1)
	remainder := queryLen % (seedLen + 1)
	// Estimate the number of gaps. This is a simplification and might need adjustment.
	// For example, assume one gap per seed for simplicity.
	estimatedGaps := seeds

	// Calculate the penalty from opening and extending gaps
	gapPenalty := estimatedGaps*config.OpenGapPenalty + (seedLen-1)*estimatedGaps*config.GapPenalty

	// Adjust the score calculations to include the gap penalty
	if seedLen*config.MaxMatch-gapPenalty >= currBestScore && perfectScore-((queryLen-seedLen)*config.MinMatch)-gapPenalty >= currBestScore {
		return true
	} else if seedLen*seeds*config.MaxMatch+seeds*config.LeastSevereMismatch-gapPenalty >= currBestScore &&
		perfectScore-remainder*config.MinMatch+seeds*config.LeastSevereMatchMismatchChange-gapPenalty >= currBestScore {
		return true
	} else if seedLen*seeds*config.MaxMatch+remainder*config.MaxMatch+(seeds+1)*config.LeastSevereMismatch-gapPenalty >= currBestScore && perfectScore+(seeds+1)*config.LeastSevereMatchMismatchChange-gapPenalty >= currBestScore {
		return true
	} else {
		return false
	}
}

func ChromAndPosToNumber(chrom int, start int) uint64 {
	return (uint64(chrom) << 32) | uint64(start)
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
	chromIdx := (code & leftMast) >> 32
	pos := code & rightMask
	return int64(chromIdx), int64(pos)
}

func SortSeedDevByTotalLen(seeds []*SeedDev) {
	sort.Slice(seeds, func(i, j int) bool { return seeds[i].TotalLength > seeds[j].TotalLength })
}
