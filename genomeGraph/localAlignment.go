package genomeGraph

import (
	"log"
	"math"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

// MatrixMemoryPool represents a scoring matrix used in local alignment, specifically for sync.Pool
type MatrixMemoryPool struct {
	i        int
	j        int
	routeIdx int

	matrix [][]int64
	trace  [][]byte

	Seq   []dna.Base
	Path  []uint32
	Route []cigar.ByteCigar

	QueryStart  int
	QueryEnd    int
	TargetStart int
	TargetEnd   int

	CurrMax   int64
	CurrScore int64
}

type ScoreCard struct {
	Curr *SeedDev
	Seq  []dna.Base
	Tail *SeedDev

	TargetStart int
	TargetEnd   int
	QueryStart  int
	QueryEnd    int

	PerfectScore int64

	CurrScore int64
	SeedScore int64

	LeftScore  int64
	RightScore int64
	LeftPath   []uint32
	RightPath  []uint32

	LeftAlignment  []cigar.ByteCigar
	RightAlignment []cigar.ByteCigar
}

// MatrixPoolMemory returns a sync.Pool that provides a pool of reusable MatrixScore objects to reduce memory allocation overhead
func MatrixPoolMemory(size int) *sync.Pool {
	m := make([][]int64, size)
	t := make([][]byte, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		t[idx] = make([]byte, size)
	}
	return &sync.Pool{
		New: func() interface{} {
			memory := MatrixMemoryPool{
				i:           0,
				j:           0,
				routeIdx:    0,
				matrix:      m,
				trace:       t,
				Seq:         make([]dna.Base, 0, 150),
				Path:        make([]uint32, 0, 10),
				Route:       make([]cigar.ByteCigar, 0, 3),
				QueryStart:  0,
				TargetStart: 0,
				TargetEnd:   0,
				QueryEnd:    0,
				CurrScore:   math.MinInt64,
				CurrMax:     0,
			}
			return &memory
		},
	}
}

func SmithWaterman(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]byte) (int64, []cigar.ByteCigar, int, int, int, int) {
	//check if size of alpha is larger than m
	var currMax int64
	var maxI int
	var maxJ int
	var i, j, routeIdx int
	//setting up the first rows and columns
	//seting up the rest of the matrix
	for i := 0; i < len(alpha)+1; i++ {
		m[i][0] = 0
	}
	for j := 0; j < len(beta)+1; j++ {
		m[0][j] = 0
	}
	for i = 1; i < len(alpha)+1; i++ {
		for j = 1; j < len(beta)+1; j++ {
			m[i][j], trace[i][j] = cigar.TraceMatrixExtension(m[i-1][j-1], m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			if m[i][j] > currMax {
				currMax = m[i][j]
				maxI = i
				maxJ = j
			}
			if m[i][j] < 0 {
				m[i][j] = 0
			}
		}
	}
	var minI, minJ int
	var curr cigar.ByteCigar
	var route []cigar.ByteCigar
	route = append(route, cigar.ByteCigar{RunLen: 0, Op: trace[maxI][maxJ]})
	for i, j, routeIdx = maxI, maxJ, 0; m[i][j] > 0; {
		if route[routeIdx].RunLen == 0 {
			route[routeIdx].RunLen = 1
			route[routeIdx].Op = trace[i][j]
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLen += 1
		} else {
			curr = cigar.ByteCigar{RunLen: 0, Op: trace[i][j]}
			route = append(route, curr)
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
			log.Fatalf("Error: unexpected traceback")
		}
		minI = i
		minJ = j
	}
	cigar.ReverseBytesCigar(route)
	return m[maxI][maxJ], route, minI, maxI, minJ, maxJ
}

func LeftLocal(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]byte) (int64, []cigar.ByteCigar, int, int, int, int) {
	//check if size of alpha is larger than m
	var i, j, routeIdx int
	initialZeroMatrix(m, len(alpha), len(beta))
	for i = 1; i < len(alpha)+1; i++ {
		for j = 1; j < len(beta)+1; j++ {
			m[i][j], trace[i][j] = cigar.TraceMatrixExtension(m[i-1][j-1], m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			if m[i][j] < 0 {
				m[i][j] = 0
			}
		}
	}
	var minI, minJ = len(alpha), len(beta)
	var route []cigar.ByteCigar
	//traceback starts in top corner
	for i, j, routeIdx = len(alpha), len(beta), 0; m[i][j] > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(route) == 0 {
			curr := cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = append(route, curr)
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLen += 1
		} else {
			curr := cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = append(route, curr)
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
	cigar.ReverseBytesCigar(route)
	return m[len(alpha)][len(beta)], route, minI, len(alpha), minJ, len(beta)
}

func RightLocal(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]byte) (int64, []cigar.ByteCigar, int, int, int, int) {
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
				m[i][j], trace[i][j] = cigar.TraceMatrixExtension(m[i-1][j-1], m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			}
			if m[i][j] > currMax {
				currMax = m[i][j]
				maxI = i
				maxJ = j
			}
		}
	}
	var route []cigar.ByteCigar
	//traceback starts in top corner
	for i, j, routeIdx = maxI, maxJ, 0; i > 0 || j > 0; {
		//if route[routeIdx].RunLength == 0 {
		if len(route) == 0 {
			curr := cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = append(route, curr)
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLen += 1
		} else {
			curr := cigar.ByteCigar{RunLen: 1, Op: trace[i][j]}
			route = append(route, curr)
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
	cigar.ReverseBytesCigar(route)
	return m[maxI][maxJ], route, 0, maxI, 0, maxJ
}

func tripleMaxTrace(prev int64, a int64, b int64, c int64) (int64, rune) {
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

func initialZeroMatrix(m [][]int64, alphaLen int, betaLen int) {
	for i := 0; i < alphaLen+1; i++ {
		m[i][0] = 0
	}
	for j := 0; j < betaLen+1; j++ {
		m[0][j] = 0
	}
}
