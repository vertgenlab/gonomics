package qDna

import (
	"log"
	//"math"
	"github.com/vertgenlab/gonomics/align"
)

//Modified from align package (linear.go) takes a sequence of qbases instead of normal bases
func ConstGap(alpha []*QBase, beta []*QBase, scores [][]float64, gapPen float64) (float64, []align.Cigar) {
	m := make([][]float64, len(alpha)+1)
	trace := make([][]align.ColType, len(alpha)+1)
	for idx := range m {
		m[idx] = make([]float64, len(beta)+1)
		trace[idx] = make([]align.ColType, len(beta)+1)
	}

	var i, j, routeIdx int
	for i = range m {
		for j = range m[0] {
			if i == 0 && j == 0 {
				m[i][j] = 0
			} else if i == 0 {
				m[i][j] = m[i][j-1] + gapPen
				trace[i][j] = 1
			} else if j == 0 {
				m[i][j] = m[i-1][j] + gapPen
				trace[i][j] = 2
			} else {
				m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1]+QDnaScore(alpha[i-1], beta[j-1], scores), m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			}
		}
	}

	route := make([]align.Cigar, 1)
	for i, j, routeIdx = len(trace)-1, len(trace[0])-1, 0; i > 0 || j > 0; {
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
		case 1:
			j -= 1
		case 2:
			i -= 1
		default:
			log.Fatalf("Error: unexpected traceback")
		}
	}
	reverseCigar(route)
	return m[len(m)-1][len(m[0])-1], route
}

func AffineGap(alpha []*QBase, beta []*QBase, scores [][]float64, gapOpen float64, gapExtend float64) (float64, []align.Cigar) {
	m, trace := initAffineScoringAndTrace(len(alpha), len(beta))
	for i := range m[0] {
		for j := range m[0][0] {
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
}

func affineTrace(m [][][]float64, trace [][][]align.ColType) (float64, []align.Cigar) {
	route := make([]align.Cigar, 1)
	lastI := len(m[0]) - 1
	lastJ := len(m[0][0]) - 1
	maxScore, k := tripleMaxTrace(m[0][lastI][lastJ], m[1][lastI][lastJ], m[2][lastI][lastJ])
	for i, j, routeIdx := lastI, lastJ, 0; i > 0 || j > 0; {
		if route[routeIdx].RunLength == 0 {
			route[routeIdx].RunLength = 1
			route[routeIdx].Op = k
		} else if route[routeIdx].Op == k {
			route[routeIdx].RunLength += 1
		} else {
			route = append(route, align.Cigar{RunLength: 1, Op: k})
			routeIdx++
		}
		switch k {
		case align.ColM:
			k = trace[k][i][j]
			i--
			j--
		case align.ColI:
			k = trace[k][i][j]
			j--
		case align.ColD:
			k = trace[k][i][j]
			i--
		default:
			log.Fatalf("Error: unexpected traceback")
		}
	}
	reverseCigar(route)
	return maxScore, route
}

func initAffineScoringAndTrace(firstSeqLen int, secondSeqLen int) ([][][]float64, [][][]align.ColType) {
	m := make([][][]float64, 3)
	trace := make([][][]align.ColType, 3)
	for k := range m {
		m[k] = make([][]float64, firstSeqLen+1)
		trace[k] = make([][]align.ColType, firstSeqLen+1)
		for i := range m[0] {
			m[k][i] = make([]float64, secondSeqLen+1)
			trace[k][i] = make([]align.ColType, secondSeqLen+1)
		}
	}
	return m, trace
}

func tripleMaxTrace(a float64, b float64, c float64) (float64, align.ColType) {
	if a >= b && a >= c {
		return a, align.ColM
	} else if b >= c {
		return b, align.ColI
	} else {
		return c, align.ColD
	}
}

//var veryNegNum float64 = float64(math.MinInt64 / 2)
