package qDna

import (
//"log"

)

func UngappedAlign(alpha []*QBase, beta []*QBase, scoreMatrix [][]float64) (float64, int64, int64, int64, int64) {
	var score, bestScore float64 = 0, 0
	var bestStartI, bestEndI, bestStartJ, bestEndJ int64 = 0, 0, 0, 0
	var i, j, startI, endI, startJ, endJ int64 = 0, 0, 0, 0, 0, 0
	for j = int64(len(beta) - 1); j > 0; j-- {
		score, startI, endI, startJ, endJ = ungappedRegionScore(alpha, 0, beta, j, scoreMatrix)
		if score > bestScore {
			bestScore = score
			bestStartI, bestEndI, bestStartJ, bestEndJ = startI, endI, startJ, endJ
		}
	}
	for i = 0; i < int64(len(alpha)); i++ {
		score, startI, endI, startJ, endJ = ungappedRegionScore(alpha, i, beta, 0, scoreMatrix)
		if score > bestScore {
			bestScore = score
			bestStartI, bestEndI, bestStartJ, bestEndJ = startI, endI, startJ, endJ
		}
	}
	return bestScore, bestStartI, bestEndI, bestStartJ, bestEndJ
}

func ungappedRegionScore(alpha []*QBase, alphaStart int64, beta []*QBase, betaStart int64, scores [][]float64) (float64, int64, int64, int64, int64) {
	var answer, currMax float64 = 0, 0
	var bestStartI, bestStartJ int64 = 0, 0
	var i, j, startI, endI, startJ, endJ int64 = 0, 0, 0, 0, 0, 0

	for i, j = alphaStart, betaStart; i < int64(len(alpha)) && j < int64(len(beta)); i, j = i+1, j+1 {
		answer += QDnaScore(alpha[i], beta[j], scores)
		if answer <= 0 {
			answer = 0
			startI = i
			startJ = j
		}
		if answer > currMax {
			currMax = answer
			bestStartI = startI
			bestStartJ = startJ
			endI = i
			endJ = j
		}

	}
	return currMax, bestStartI, endI, bestStartJ, endJ
}
