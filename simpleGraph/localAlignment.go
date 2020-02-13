package simpleGraph

import (
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	//"fmt"
)

var HumanChimpTwoScoreMatrixNoGap = [][]int64{
	{90, -330, -236, -356},
	{-330, 100, -318, -236},
	{-236, -318, 100, -330},
	{-356, -236, -330, 90},
}

func reverseCigar(alpha []align.Cigar) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
}

func reverseCigarPointer(alpha []*cigar.Cigar) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
}

func swMatrixSetup(size int64) ([][]int64, [][]rune) {
	m := make([][]int64, size)
	trace := make([][]rune, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		trace[idx] = make([]rune, size)
	}
	return m, trace
}

func SwMatrixSetup(size int64) ([][]int64, [][]rune) {
	m := make([][]int64, size)
	trace := make([][]rune, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		trace[idx] = make([]rune, size)
	}
	return m, trace
}

func initialZeroMatrix(m [][]int64, alphaLen int, betaLen int) {
	for i := 0; i < alphaLen+1; i++ {
		m[i][0] = 0
	}
	for j := 0; j < betaLen+1; j++ {
		m[0][j] = 0
	}
}

func initialTraceMatrix(trace [][]rune, alphaLen int, betaLen int) {
	for i := 1; i < alphaLen+1; i++ {
		trace[i][0] = 'D'
	}
	for j := 1; j < betaLen+1; j++ {
		trace[0][j] = 'I'
	}
}

/*
func initializeLeftGlobal(m [][]int64, trace [][]rune, alphaLen int, betaLen int, gapPen int64) {
	for i:=0, i < alphaLen+1; i++ {
                for j:=0; j<betaLen+1; j++ {
                        if i == 0 && j == 0 {
                                m[i][j] = 0
                        } else if i == 0 {
                                m[i][j] = m[i][j-1] + gapPen
                                trace[i][j] = 'I'
                        } else if j == 0 {
                                m[i][j] = m[i-1][j] + gapPen
                                trace[i][j] = 'D'
                        } else {
                                m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
                        }
                }
        }
}*/

func SmithWaterman(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]rune) (int64, []*cigar.Cigar, int64, int64, int64, int64) {
	//check if size of alpha is larger than m
	var currMax int64
	var maxI int64
	var maxJ int64
	var i, j, routeIdx int64
	//setting up the first rows and columns
	//seting up the rest of the matrix
	initialZeroMatrix(m, len(alpha), len(beta))
	for i = 1; i < int64(len(alpha)+1); i++ {
		for j = 1; j < int64(len(beta)+1); j++ {
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
	var minI, minJ int64
	var curr cigar.Cigar
	var route []*cigar.Cigar
	route = append(route, &cigar.Cigar{RunLength: 0, Op: trace[maxI][maxJ]})
	for i, j, routeIdx = maxI, maxJ, 0; m[i][j] > 0; {
		if route[routeIdx].RunLength == 0 {
			route[routeIdx].RunLength = 1
			route[routeIdx].Op = trace[i][j]
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLength += 1
		} else {
			curr = cigar.Cigar{RunLength: 0, Op: trace[i][j]}
			route = append(route, &curr)
			routeIdx++
		}
		switch trace[i][j] {
		case 'M':
			i, j = i-1, j-1
			//refStart = refStart + 1
		case 'I':
			j -= 1
		case 'D':
			i -= 1
			//refStart = refStart + 1
		default:
			log.Fatalf("Error: unexpected traceback")
		}
		minI = i
		minJ = j
	}
	reverseCigarPointer(route)
	return m[maxI][maxJ], route, minI, maxI, minJ, maxJ
}

func tripleMaxTrace(a int64, b int64, c int64) (int64, rune) {
	if a >= b && a >= c {
		return a, 'M'
	} else if b >= c {
		return b, 'I'
	} else {
		return c, 'D'
	}
}

func LeftLocal(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]rune) (int64, []*cigar.Cigar, int, int, int, int) {
	//check if size of alpha is larger than m
	//var currMax int64
	//var maxI int64
	//var maxJ int64
	var i, j, routeIdx int
	initialZeroMatrix(m, len(alpha), len(beta))
	for i = 1; i < len(alpha)+1; i++ {
		for j = 1; j < len(beta)+1; j++ {
			m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			/*if m[i][j] > currMax {
				currMax = m[i][j]
				maxI = i
				maxJ = j
			}*/
			if m[i][j] < 0 {
				m[i][j] = 0
			}
		}
	}
	var minI, minJ = len(alpha), len(beta)

	var route []*cigar.Cigar
	//route = append(route, &cigar.Cigar{RunLength: 0, Op: trace[len(alpha)][len(beta)]})
	//traceback starts in top corner
	//fmt.Println("a little more!")
	for i, j, routeIdx = len(alpha), len(beta), 0; m[i][j] > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(route) == 0 {
			curr := cigar.Cigar{RunLength: 1, Op: trace[i][j]}
			route = append(route, &curr)
			//route[routeIdx].RunLength = 1
			//route[routeIdx].Op = trace[i][j]
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLength += 1
		} else {
			curr := cigar.Cigar{RunLength: 1, Op: trace[i][j]}
			route = append(route, &curr)
			routeIdx++
		}
		switch trace[i][j] {
		case 'M':
			i, j = i-1, j-1
		case 'I':
			j -= 1
		case 'D':
			i -= 1
		default:
			log.Fatalf("Error: unexpected traceback")
		}
		minI = i
		minJ = j
		//log.Printf("cigar in the making=%s\n", cigar.ToString(route))
	}
	//TODO: double check if this is tracing back in the correct directions
	reverseCigarPointer(route)
	return m[len(alpha)][len(beta)], route, minI, len(alpha), minJ, len(beta)
}

func RightLocal(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]rune) (int64, []*cigar.Cigar, int, int, int, int) {
	//check if size of alpha is larger than m
	var currMax int64
	var maxI int
	var maxJ int
	var i, j, routeIdx int

	//setting up the first rows and columns
	//seting up the rest of the matrix

	for i = 0; i < len(alpha)+1; i++ {
		for j = 0; j < len(beta)+1; j++ {
			if i == 0 && j == 0 {
				m[i][j] = 0
			} else if i == 0 {
				m[i][j] = m[i][j-1] + gapPen
				trace[i][j] = 'I'
			} else if j == 0 {
				m[i][j] = m[i-1][j] + gapPen
				trace[i][j] = 'D'
			} else {
				m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			}
			if m[i][j] > currMax {
				currMax = m[i][j]
				maxI = i
				maxJ = j
			}
		}
	}

	//var minI, minJ int64 = 1, 1
	//var refStart int
	var route []*cigar.Cigar
	//route = append(route, &cigar.Cigar{RunLength: 0, Op: trace[maxI][maxJ]})
	//traceback starts in top corner
	for i, j, routeIdx = maxI, maxJ, 0; i > 0 || j > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(route) == 0 {
			curr := cigar.Cigar{RunLength: 1, Op: trace[i][j]}
			route = append(route, &curr)
			//route[routeIdx].RunLength = 1
			//route[routeIdx].Op = trace[i][j]
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLength += 1
		} else {
			curr := cigar.Cigar{RunLength: 1, Op: trace[i][j]}
			route = append(route, &curr)
			routeIdx++
		}
		switch trace[i][j] {
		case 'M':
			i, j = i-1, j-1
		//	refStart = refStart + 1
		case 'I':
			j -= 1
		case 'D':
			i -= 1
		//	refStart = refStart + 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", trace[i][j])
		}
	}
	reverseCigarPointer(route)

	return m[maxI][maxJ], route, 0, maxI, 0, maxJ
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
