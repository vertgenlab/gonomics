package genomeGraph

import (
	"log"
	"math"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

func reverseCigarPointer(alpha []cigar.Cigar) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
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

func initialTraceMatrix(trace [][]rune, alphaLen int, betaLen int) {
	for i := 1; i < alphaLen+1; i++ {
		trace[i][0] = 'D'
	}
	for j := 1; j < betaLen+1; j++ {
		trace[0][j] = 'I'
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

type MatrixScore struct {
	i        int
	j        int
	routeIdx int

	matrix [][]int64
	trace  [][]byte

	Seq         []dna.Base
	Path        []uint32
	route       []cigar.ByteCigar
	queryStart  int
	queryEnd    int
	targetStart int
	targetEnd   int

	currMax   int64
	currScore int64
}

type GraphSettings struct {
	scores     [][]int64
	gapPenalty int64
	tileSize   int
	stepSize   int
	extension  int
}

func MatrixPoolMemory() *sync.Pool {
	m := make([][]int64, defaultMatrixSize)
	t := make([][]byte, defaultMatrixSize)
	for idx := range m {
		m[idx] = make([]int64, defaultMatrixSize)
		t[idx] = make([]byte, defaultMatrixSize)
	}
	return &sync.Pool{
		New: func() interface{} {
			memory := MatrixScore{
				i:           0,
				j:           0,
				routeIdx:    0,
				matrix:      m,
				trace:       t,
				Seq:         make([]dna.Base, 0, 150),
				Path:        make([]uint32, 0, 10),
				route:       make([]cigar.ByteCigar, 0, 3),
				queryStart:  0,
				targetStart: 0,
				targetEnd:   0,
				queryEnd:    0,
				currScore:   math.MinInt64,
				currMax:     0,
			}
			return &memory
		},
	}
}
