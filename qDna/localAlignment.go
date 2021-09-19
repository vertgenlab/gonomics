package qDna

import (
	//"fmt"
	"log"
	//"math"
	"github.com/vertgenlab/gonomics/align"
	//"github.com/vertgenlab/gonomics/common"
)

func SmithWaterman(alpha []*QBase, beta []*QBase, scores [][]float64, gapPen float64) (float64, []align.Cigar, int64, int64, int64, int64) {

	m := make([][]float64, len(alpha)+1)
	trace := make([][]align.ColType, len(alpha)+1)
	for idx := range m {
		m[idx] = make([]float64, len(beta)+1)
		trace[idx] = make([]align.ColType, len(beta)+1)
	}
	var currMax float64
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
			m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1]+QDnaScore(alpha[i-1], beta[j-1], scores), m[i][j-1]+gapPen, m[i-1][j]+gapPen)
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
	//for x := 0; x < len(m); x++ {
	//	for y := 0; y < len(m[0]); y ++ {
	//		fmt.Print(m[x][y], ", ")
	//	}
	//	fmt.Println("")
	//}
	//fmt.Println(align.View(mostLikelySeq(alpha), mostLikelySeq(beta), maxI))
	return m[maxI][maxJ], route, minI, maxI, minJ, maxJ
}

/*
func AffineGapSW(alpha []*QBase, beta []*QBase, scores [][]float64, gapOpen float64, gapExtend float64) (float64, []align.Cigar) {
	m, trace := initAffineScoringAndTrace(len(alpha), len(beta))
	for i, _ := range m[0] {
		for j, _ := range m[0][0] {
			if i == 0 && j == 0 {
				m[0][i][j] = 0
				m[1][i][j] = gapOpen
				m[2][i][j] = gapOpen
			} else if i == 0 {
				//m[0][i][j] = veryNegNum
				m[0][i][j] = 0
				m[1][i][j] = gapExtend + m[1][i][j-1]
				trace[1][i][j] = align.ColI
				//m[2][i][j] = veryNegNum
				m[2][i][j] = 0
			} else if j == 0 {
				//m[0][i][j] = veryNegNum
				//m[1][i][j] = veryNegNum
				m[0][i][j] = 0
				m[1][i][j] = 0
				m[2][i][j] = gapExtend + m[2][i-1][j]
				trace[2][i][j] = align.ColD
			} else {
				m[0][i][j], trace[0][i][j] = tripleMaxTrace(QDnaScore(alpha[i-1], beta[j-1], scores)+m[0][i-1][j-1], QDnaScore(alpha[i-1], beta[j-1], scores)+m[1][i-1][j-1], QDnaScore(alpha[i-1], beta[j-1], scores)+m[2][i-1][j-1])
				m[1][i][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend+m[0][i][j-1], gapExtend+m[1][i][j-1], gapOpen+gapExtend+m[2][i][j-1])
				m[2][i][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend+m[0][i-1][j], gapOpen+gapExtend+m[1][i-1][j], gapExtend+m[2][i-1][j])
			}
		}
	}
	maxScore, route := affineTrace(m, trace)

	return maxScore, route
}*/
