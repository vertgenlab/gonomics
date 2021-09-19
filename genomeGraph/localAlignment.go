package genomeGraph

import (
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

var HumanChimpTwoScoreMatrixNoGap = [][]int{
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

func swMatrixSetup(size int) ([][]int, [][]rune) {
	m := make([][]int, size)
	trace := make([][]rune, size)
	for idx := range m {
		m[idx] = make([]int, size)
		trace[idx] = make([]rune, size)
	}
	return m, trace
}

func initialZeroMatrix(m [][]int, alphaLen int, betaLen int) {
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

func SmithWaterman(alpha []dna.Base, beta []dna.Base, scores [][]int, gapPen int, m [][]int, trace [][]rune) (int, []*cigar.Cigar, int, int, int, int) {
	//check if size of alpha is larger than m
	var currMax int
	var maxI int
	var maxJ int
	var i, j, routeIdx int
	//setting up the first rows and columns
	//seting up the rest of the matrix
	initialZeroMatrix(m, len(alpha), len(beta))
	for i = 1; i < int(len(alpha)+1); i++ {
		for j = 1; j < int(len(beta)+1); j++ {
			m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1], m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			if m[i][j] > currMax {
				currMax = m[i][j]
				maxI = int(i)
				maxJ = int(j)
			}
			if m[i][j] < 0 {
				m[i][j] = 0
			}
		}
	}
	var minI, minJ int
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
		case '=':
			i, j = i-1, j-1
			//refStart = refStart + 1
		case 'X':
			i, j = i-1, j-1
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

func tripleMaxTrace(prev int, a int, b int, c int) (int, rune) {
	if a >= b && a >= c {
		if a > prev {
			return a, '='
		} else {
			return a, 'X'
		}
	} else if b >= c {
		return b, 'I'
	} else {
		return c, 'D'
	}
}

func LeftLocal(alpha []dna.Base, beta []dna.Base, scores [][]int, gapPen int, m [][]int, trace [][]rune) (int, []*cigar.Cigar, int, int, int, int) {
	//check if size of alpha is larger than m
	var i, j, routeIdx int
	initialZeroMatrix(m, len(alpha), len(beta))
	for i = 1; i < len(alpha)+1; i++ {
		for j = 1; j < len(beta)+1; j++ {
			m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1], m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			if m[i][j] < 0 {
				m[i][j] = 0
			}
		}
	}
	var minI, minJ = len(alpha), len(beta)
	var route []*cigar.Cigar
	//traceback starts in top corner
	for i, j, routeIdx = len(alpha), len(beta), 0; m[i][j] > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(route) == 0 {
			curr := cigar.Cigar{RunLength: 1, Op: trace[i][j]}
			route = append(route, &curr)
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLength += 1
		} else {
			curr := cigar.Cigar{RunLength: 1, Op: trace[i][j]}
			route = append(route, &curr)
			routeIdx++
		}
		switch trace[i][j] {
		case '=':
			i, j = i-1, j-1
		case 'X':
			i, j = i-1, j-1
		case 'I':
			j -= 1
		case 'D':
			i -= 1
		default:
			log.Fatalf("Error: unexpected traceback %c\n", trace[i][j])
		}
		minI = i
		minJ = j
	}
	//TODO: double check if this is tracing back in the correct directions
	reverseCigarPointer(route)
	return m[len(alpha)][len(beta)], route, minI, len(alpha), minJ, len(beta)
}

func RightLocal(alpha []dna.Base, beta []dna.Base, scores [][]int, gapPen int, m [][]int, trace [][]rune) (int, []*cigar.Cigar, int, int, int, int) {
	//check if size of alpha is larger than m
	var currMax int
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
				m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1], m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			}
			if m[i][j] > currMax {
				currMax = m[i][j]
				maxI = i
				maxJ = j
			}
		}
	}
	var route []*cigar.Cigar
	//traceback starts in top corner
	for i, j, routeIdx = maxI, maxJ, 0; i > 0 || j > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(route) == 0 {
			curr := cigar.Cigar{RunLength: 1, Op: trace[i][j]}
			route = append(route, &curr)
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLength += 1
		} else {
			curr := cigar.Cigar{RunLength: 1, Op: trace[i][j]}
			route = append(route, &curr)
			routeIdx++
		}
		switch trace[i][j] {
		case '=':
			i, j = i-1, j-1
		case 'X':
			i, j = i-1, j-1
		case 'I':
			j -= 1
		case 'D':
			i -= 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", trace[i][j])
		}
	}
	reverseCigarPointer(route)
	return m[maxI][maxJ], route, 0, maxI, 0, maxJ
}
