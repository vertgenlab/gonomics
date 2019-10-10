package simpleGraph

import (
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

var HumanChimpTwoScoreMatrixNoGap = [][]int64{
	{90, -330, -236, -356},
	{-330, 100, -318, -236},
	{-236, -318, 100, -330},
	{-356, -236, -330, 90},
}

func SmithWaterman(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64) (int64, []align.Cigar, int64, int64, int64, int64) {
	m := make([][]int64, len(alpha)+1)
	trace := make([][]align.ColType, len(alpha)+1)
	for idx := range m {
		m[idx] = make([]int64, len(beta)+1)
		trace[idx] = make([]align.ColType, len(beta)+1)
	}
	var currMax int64
	var maxI int64
	var maxJ int64
	var i, j, routeIdx int64

	//setting up the first rows and columns
	for i = 0; i < int64(len(m)); i++ {
		m[i][0] = 0
	}
	for j = 0; j < int64(len(m[0])); j++ {
		m[0][j] = 0
	}
	//seting up the rest of the matrix
	//var diagScore, upScore, leftScore float64
	for i = 1; i < int64(len(m)); i++ {
		for j = 1; j < int64(len(m[0])); j++ {
			m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			if m[i][j] > currMax {
				currMax = m[i][j]
				maxI = int64(i)
				maxJ = int64(j)
			}
			if m[i][j] < 0 {
				m[i][j] = 0
			}
		}
	}
	route := make([]align.Cigar, 1)
	var minI, minJ int64
	var refStart int
	for i, j, routeIdx = maxI, maxJ, 0; m[i][j] > 0; {
		if route[routeIdx].RunLength == 0 {
			route[routeIdx].RunLength = 1
			route[routeIdx].Op = trace[i][j]
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLength += 1
		} else {
			route = append(route, align.Cigar{RunLength: 1, Op: trace[i][j]})
			routeIdx++
		}
		switch trace[i][j] {
		case 0:
			i, j = i-1, j-1
			refStart = refStart + 1
		case 1:
			j -= 1
		case 2:
			i -= 1
			refStart = refStart + 1
		default:
			log.Fatalf("Error: unexpected traceback")
		}
		minI = i
		minJ = j
	}
	reverseCigar(route)
	//fmt.Println(align.LocalView(mostLikelySeq(alpha), mostLikelySeq(beta), route, maxI))
	return m[maxI][maxJ], route, minI, maxI, minJ, maxJ
}

func reverseCigar(alpha []align.Cigar) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
}

func tripleMaxTrace(a int64, b int64, c int64) (int64, align.ColType) {
	if a >= b && a >= c {
		return a, align.ColM
	} else if b >= c {
		return b, align.ColI
	} else {
		return c, align.ColD
	}
}

/*
func UngappedAlign(alpha []dna.Base, beta []dna.Base, scoreMatrix [][]int64) (int64, int64, int64, int64, int64) {
	var score, bestScore int64 = 0, 0
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

func ungappedRegionScore(alpha []dna.Base, alphaStart int64, beta []dna.Base, betaStart int64, scores [][]int64) (int64, int64, int64, int64, int64) {
	var answer, currMax int64 = 0, 0
	var bestStartI, bestStartJ int64 = 0, 0
	var i, j, startI, endI, startJ, endJ int64 = 0, 0, 0, 0, 0, 0

	for i, j = alphaStart, betaStart; i < int64(len(alpha)) && j < int64(len(beta)); i, j = i+1, j+1 {
		//answer += QDnaScore(alpha[i], beta[j], scores)
		answer += scores[alpha[i]][beta[j]]
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
}*/
