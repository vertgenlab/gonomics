package genomeGraph

import (
	"log"
	"sort"

	"github.com/vertgenlab/gonomics/dna"
)

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

func MismatchStats(scoreMatrix [][]int64) (int64, int64, int64, int64) {
	var maxMatch int64 = 0
	var minMatch int64
	var leastSevereMismatch int64 = scoreMatrix[0][1]
	var i, j int
	for i = 0; i < len(scoreMatrix); i++ {
		for j = 0; j < len(scoreMatrix[i]); j++ {
			if scoreMatrix[i][j] > maxMatch {
				minMatch = maxMatch
				maxMatch = scoreMatrix[i][j]
			} else {
				if scoreMatrix[i][j] < 0 && leastSevereMismatch < scoreMatrix[i][j] {
					leastSevereMismatch = scoreMatrix[i][j]
				}
			}
		}
	}
	var leastSevereMatchMismatchChange int64 = leastSevereMismatch - maxMatch
	return maxMatch, minMatch, leastSevereMismatch, leastSevereMatchMismatchChange
}

type ScoreMatrixHelper struct {
	Matrix                         [][]int64
	MaxMatch                       int64
	MinMatch                       int64
	LeastSevereMismatch            int64
	LeastSevereMatchMismatchChange int64
}

func getScoreMatrixHelp(scoreMatrix [][]int64) *ScoreMatrixHelper {
	help := ScoreMatrixHelper{Matrix: scoreMatrix}
	help.MaxMatch, help.MinMatch, help.LeastSevereMismatch, help.LeastSevereMatchMismatchChange = MismatchStats(scoreMatrix)
	return &help
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

func SortSeedDevByTotalLen(seeds []*Seed) {
	sort.Slice(seeds, func(i, j int) bool { return seeds[i].TotalLength > seeds[j].TotalLength })
}
