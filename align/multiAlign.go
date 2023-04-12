package align

import (
	"fmt"
	"math"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

func fastaListToIndividualGroups(records []fasta.Fasta) [][]fasta.Fasta {
	answer := make([][]fasta.Fasta, len(records))
	for i := range answer {
		answer[i] = make([]fasta.Fasta, 1)
		answer[i][0] = records[i]
	}
	return answer
}

func mergeFastaGroups(groups [][]fasta.Fasta, x int, y int, route []Cigar) [][]fasta.Fasta {
	groups[x] = mergeMultipleAlignments(groups[x], groups[y], route)
	groups[y] = groups[len(groups)-1]
	groups = groups[:len(groups)-1]
	return groups
}

func nearestGroups(groups [][]fasta.Fasta, scoreMatrix [][]int64, gapOpen int64, gapExtend int64) (bestX int, bestY int, bestScore int64, bestRoute []Cigar) {
	var route []Cigar
	var score int64
	bestScore = math.MinInt64
	for x := 0; x < len(groups)-1; x++ {
		for y := x + 1; y < len(groups); y++ {
			score, route = multipleAffineGap(groups[x], groups[y], scoreMatrix, gapOpen, gapExtend)
			if score > bestScore {
				bestX, bestY, bestScore, bestRoute = x, y, score, route
			}
		}
	}
	return bestX, bestY, bestScore, bestRoute
}

func nearestGroupsChunk(groups [][]fasta.Fasta, scoreMatrix [][]int64, gapOpen int64, gapExtend int64, chunkSize int) (bestX int, bestY int, bestScore int64, bestRoute []Cigar) {
	var route []Cigar
	var score int64
	bestScore = math.MinInt64
	for x := 0; x < len(groups)-1; x++ {
		for y := x + 1; y < len(groups); y++ {
			score, route = multipleAffineGapChunk(groups[x], groups[y], scoreMatrix, gapOpen, gapExtend, int64(chunkSize))
			if score > bestScore {
				bestX, bestY, bestScore, bestRoute = x, y, score, route
			}
		}
	}
	return bestX, bestY, bestScore, bestRoute
}

func AllSeqAffine(records []fasta.Fasta, scoreMatrix [][]int64, gapOpen int64, gapExtend int64) []fasta.Fasta {
	groups := fastaListToIndividualGroups(records)
	for len(groups) > 1 {
		x, y, _, route := nearestGroups(groups, scoreMatrix, gapOpen, gapExtend)
		groups = mergeFastaGroups(groups, x, y, route)
	}
	return groups[0]
}

//align sequences
func AllSeqAffineChunk(records []fasta.Fasta, scoreMatrix [][]int64, gapOpen int64, gapExtend int64, chunkSize int) []fasta.Fasta {
	groups := fastaListToIndividualGroups(records)
	for len(groups) > 1 {
		x, y, score, route := nearestGroupsChunk(groups, scoreMatrix, gapOpen, gapExtend, chunkSize)
		fmt.Printf("x=%d ; y=%d ; score=%d ; cigar=%v ; len(groups)=%d\n", x, y, score, route, len(groups))
		groups = mergeFastaGroups(groups, x, y, route)
	}
	return groups[0]
}

//average of pairs scoring scheme where gaps are ignored
//maybe there should be a small penalty for gaps so that gaps will tend to be in the same location
func scoreColumnMatch(alpha []fasta.Fasta, beta []fasta.Fasta, alphaCol int, betaCol int, scores [][]int64) int64 {
	var sum, count int64 = 0, 0
	var alphaBase, betaBase dna.Base
	for alphaSeqIdx := range alpha {
		alphaBase = alpha[alphaSeqIdx].Seq[alphaCol]
		if alphaBase >= 5 && alphaBase <= 9 { // converts to uppercase base
			alphaBase -= 5
		}
		for betaSeqIdx := range beta {
			betaBase = beta[betaSeqIdx].Seq[betaCol]
			if betaBase >= 5 && betaBase <= 9 { // converts to uppercase base
				betaBase -= 5
			}
			if alphaBase != dna.Gap && betaBase != dna.Gap {
				sum += scores[alphaBase][betaBase]
				count++
			}
		}
	}
	return sum / count
}

func ungappedRegionColumnScore(alpha []fasta.Fasta, alphaStart int, beta []fasta.Fasta, betaStart int, length int, scores [][]int64) int64 {
	var answer int64 = 0
	for i, j := alphaStart, betaStart; i < alphaStart+length; i, j = i+1, j+1 {
		answer += scoreColumnMatch(alpha, beta, i, j, scores)
	}
	return answer
}

func mergeMultipleAlignments(alpha []fasta.Fasta, beta []fasta.Fasta, route []Cigar) []fasta.Fasta {
	answer := make([]fasta.Fasta, len(alpha)+len(beta))
	totalCols := countAlignmentColumns(route)

	for i := range answer {
		if i < len(alpha) {
			answer[i] = fasta.Fasta{Name: alpha[i].Name, Seq: make([]dna.Base, totalCols)}
		} else {
			answer[i] = fasta.Fasta{Name: beta[i-len(alpha)].Name, Seq: make([]dna.Base, totalCols)}
		}
	}

	var alphaCol, betaCol, ansCol int = 0, 0, 0
	for i := range route {
		for j := 0; j < int(route[i].RunLength); j++ {
			for k := range answer {
				if k < len(alpha) {
					if route[i].Op == ColM || route[i].Op == ColD {
						answer[k].Seq[ansCol] = alpha[k].Seq[alphaCol]
					} else {
						answer[k].Seq[ansCol] = dna.Gap
					}
				} else {
					if route[i].Op == ColM || route[i].Op == ColI {
						answer[k].Seq[ansCol] = beta[k-len(alpha)].Seq[betaCol]
					} else {
						answer[k].Seq[ansCol] = dna.Gap
					}
				}
			}
			switch route[i].Op {
			case ColM:
				alphaCol, betaCol, ansCol = alphaCol+1, betaCol+1, ansCol+1
			case ColI:
				betaCol, ansCol = betaCol+1, ansCol+1
			case ColD:
				alphaCol, ansCol = alphaCol+1, ansCol+1
			}
		}
	}
	return answer
}
