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

func swMatrixSetup(size int64) ([][]int64, [][]rune) {
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

// scores [][]int64
// gapPen int64, dynamicScore dynamicScoreKeeper
// matrix *MatrixAln,

// type dynamicScoreKeeper struct {
// 	i        int
// 	j        int
// 	routeIdx int
// 	currMax  int64
// 	route    []cigar.ByteCigar
// }

func NewMatrixPool() *sync.Pool {
	matrix, trace := MatrixSetup(defaultMatrixSize)
	return &sync.Pool{
		New: func() interface{} {
			memory := MatrixAln{
				m:     matrix,
				trace: trace,
			}
			return &memory
		},
	}
}

type SearchSettings struct {
	scoreMatrix [][]int64
	gapPenalty  int64
	tileSize    int
	extension   int
}

type SearchCache struct {
	Seq  []dna.Base
	Path []uint32

	bestScore int64
	route     []cigar.ByteCigar
	matrix    [][]int64
	trace     [][]byte

	targetStart int
	targetEnd   int
	queryStart  int
	queryEnd    int
}

func NewCache(size int) sync.Pool {
	m := make([][]int64, size)
	t := make([][]byte, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		t[idx] = make([]byte, size)
	}
	return sync.Pool{
		New: func() interface{} {
			visited := SearchCache{
				Seq:         make([]dna.Base, 0, 150),
				Path:        make([]uint32, 0, 10),
				route:       make([]cigar.ByteCigar, 0, 10),
				matrix:      m,
				trace:       t,
				targetStart: 0,
				targetEnd:   0,
				queryStart:  0,
				queryEnd:    0,
			}
			return &visited
		},
	}
}

func BreathSearchRight(n *Node, seq []dna.Base, start int, currentPath []uint32, pool *sync.Pool, settings SearchSettings) (int64, []cigar.ByteCigar, int, int, []uint32) {
	visited := pool.Get().(*SearchCache)
	defer pool.Put(visited)

	visited.Seq, visited.Path = visited.Seq[:0], visited.Path[:0]
	visited.Seq = getTargetBases(n, settings.extension, start, seq, visited.Seq, right)

	// Reuse the currentPath slice if possible
	if cap(visited.Path) >= len(currentPath) {
		visited.Path = visited.Path[:len(currentPath)]
	} else {
		visited.Path = make([]uint32, len(currentPath))
	}
	copy(visited.Path, currentPath)

	if len(seq)+len(n.Seq)-start >= settings.extension || len(n.Next) == 0 {
		visited.bestScore, visited.route, visited.targetEnd, visited.queryEnd = RightGraphAlign(visited.Seq, seq, *visited, settings)

		return visited.bestScore, visited.route, visited.targetEnd + start, visited.queryEnd, visited.Path
	} else {
		for _, i := range n.Next {
			currScore, route, targetEnd, queryEnd, path := BreathSearchRight(i.Dest, visited.Seq, 0, visited.Path, pool, settings)
			if currScore > visited.bestScore {
				visited.bestScore = currScore
				visited.route = route
				visited.targetEnd = targetEnd
				visited.queryEnd = queryEnd
				visited.Path = path
			}
		}
		cigar.ReverseBytesCigar(visited.route)
		return visited.bestScore, visited.route, visited.targetEnd + start, visited.queryEnd, visited.Path
	}
}

func BreathSearchLeft(n *Node, seq []dna.Base, position int, currentPath []uint32, pool *sync.Pool, settings SearchSettings) (int64, []cigar.ByteCigar, int, int, []uint32) {
	visited := pool.Get().(*SearchCache)
	defer pool.Put(visited)

	visited.Seq, visited.Path = visited.Seq[:0], visited.Path[:0]
	visited.Seq = getTargetBases(n, settings.extension, position, seq, visited.Seq, left)
	visited.Path = make([]uint32, len(currentPath))
	copy(visited.Path, currentPath)
	AddPath(visited.Path, n.Id)

	currSeqLen := position - len(visited.Seq) - len(seq)
	if len(seq)+position >= settings.extension || len(n.Prev) == 0 {
		visited.bestScore, visited.route, visited.targetStart, visited.queryStart = LeftGraphAlign(visited.Seq, seq, *visited, settings)

		visited.targetStart = currSeqLen + visited.targetStart
		// visited.Path = visited.Path

		return visited.bestScore, visited.route, visited.targetStart, visited.queryStart, visited.Path
	} else {
		//A very negative number
		visited.bestScore = math.MinInt64
		for _, i := range n.Prev {
			currScore, route, targetStart, queryStart, leftPath := BreathSearchLeft(i.Dest, visited.Seq, len(i.Dest.Seq), visited.Path, pool, settings)
			if currScore > visited.bestScore {
				visited.bestScore = currScore
				visited.route = route
				visited.targetStart = currSeqLen + targetStart
				visited.queryStart = queryStart
				visited.Path = leftPath
			}
		}

		cigar.ReverseBytesCigar(visited.route)
		ReversePath(visited.Path)
		return visited.bestScore, visited.route, visited.targetStart, visited.queryStart, visited.Path
	}
}
func RightGraphAlign(alpha []dna.Base, beta []dna.Base, visited SearchCache, settings SearchSettings) (int64, []cigar.ByteCigar, int, int) {
	rows, columns := len(alpha)+1, len(beta)+1
	visited.targetEnd, visited.queryEnd = 0, 0
	visited.matrix[0][0] = 0
	visited.bestScore = math.MinInt64
	var i, j int

	for j = 0; j < columns; j++ {
		visited.matrix[0][j] = int64(j) * settings.gapPenalty
		visited.trace[0][j] = cigar.Insertion
	}
	for i = 0; i < rows; i++ {
		visited.matrix[i][0] = int64(i) * settings.gapPenalty
		visited.trace[i][0] = cigar.Deletion
	}

	for i = 1; i < rows; i++ {
		for j = 1; j < columns; j++ {
			visited.matrix[i][j], visited.trace[i][j] = cigar.ByteMatrixTrace(
				visited.matrix[i-1][j-1]+settings.scoreMatrix[alpha[i-1]][beta[j-1]],
				visited.matrix[i][j-1]+settings.gapPenalty,
				visited.matrix[i-1][j]+settings.gapPenalty,
			)
			if visited.matrix[i][j] > visited.bestScore {
				visited.bestScore = visited.matrix[i][j]
				visited.targetEnd = i
				visited.queryEnd = j
			}
		}
	}
	var routeIdx int
	var route []cigar.ByteCigar
	for i, j, routeIdx = visited.targetEnd, visited.queryEnd, 0; i > 0 || j > 0; {
		if len(route) == 0 {
			route = cigar.AddCigarByte(route, cigar.ByteCigar{RunLen: 1, Op: visited.trace[i][j]})
		} else if route[routeIdx].Op == visited.trace[i][j] {
			route[routeIdx].RunLen += 1
		} else {
			route = cigar.AddCigarByte(route, cigar.ByteCigar{RunLen: 1, Op: visited.trace[i][j]})
			routeIdx++
		}
		switch visited.trace[i][j] {
		case cigar.Match:
			i, j = i-1, j-1
		case cigar.Insertion:
			j -= 1
		case cigar.Deletion:
			i -= 1
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", visited.trace[i][j])
		}
	}
	return visited.matrix[visited.targetEnd][visited.queryEnd], route, visited.targetEnd, visited.queryEnd
}

func LeftGraphAlign(alpha []dna.Base, beta []dna.Base, visited SearchCache, settings SearchSettings) (int64, []cigar.ByteCigar, int, int) {
	rows, columns := len(alpha), len(beta)
	var i, j int
	for i := 0; i < rows; i++ {
		visited.matrix[i][0] = 0
	}

	for j := 0; j < columns; j++ {
		visited.matrix[0][j] = 0
	}

	for i := 1; i < rows+1; i++ {
		for j := 1; j < columns+1; j++ {
			visited.matrix[i][j], visited.trace[i][j] = cigar.ByteMatrixTrace(
				visited.matrix[i-1][j-1]+settings.scoreMatrix[alpha[i-1]][beta[j-1]],
				visited.matrix[i][j-1]+settings.gapPenalty,
				visited.matrix[i-1][j]+settings.gapPenalty,
			)
			if visited.matrix[i][j] < 0 {
				visited.matrix[i][j] = 0
			}
		}
	}
	var routeIdx int
	var route []cigar.ByteCigar
	for i, j, routeIdx = rows, columns, 0; visited.matrix[i][j] > 0; {
		if len(route) == 0 {
			route = append(route, cigar.ByteCigar{RunLen: 1, Op: visited.trace[i][j]})
		} else if route[routeIdx].Op == visited.trace[i][j] {
			route[routeIdx].RunLen += 1
		} else {
			route = append(route, cigar.ByteCigar{RunLen: 1, Op: visited.trace[i][j]})
			routeIdx++
		}
		switch visited.trace[i][j] {
		case cigar.Match:
			i, j = i-1, j-1
		case cigar.Insertion:
			j -= 1
		case cigar.Deletion:
			i -= 1
		default:
			log.Fatalf("Error: unexpected traceback %c\n", visited.trace[i][j])
		}
	}
	return visited.matrix[rows][columns], route, i, j
}
