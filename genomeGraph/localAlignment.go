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

func initialZeroMatrix(m [][]int64, alphaLen int, betaLen int) {
	for i := 0; i < alphaLen+1; i++ {
		m[i][0] = 0
	}
	for j := 0; j < betaLen+1; j++ {
		m[0][j] = 0
	}
}

func SmithWaterman(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, m [][]int64, trace [][]rune) (int64, []cigar.Cigar, int64, int64, int64, int64) {
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
			m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1], m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
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
	var route []cigar.Cigar
	route = append(route, cigar.Cigar{RunLength: 0, Op: trace[maxI][maxJ]})
	for i, j, routeIdx = maxI, maxJ, 0; m[i][j] > 0; {
		if route[routeIdx].RunLength == 0 {
			route[routeIdx].RunLength = 1
			route[routeIdx].Op = trace[i][j]
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLength += 1
		} else {
			curr = cigar.Cigar{RunLength: 0, Op: trace[i][j]}
			route = append(route, curr)
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

func reverseCigarPointer(alpha []cigar.Cigar) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
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
