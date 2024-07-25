package genomeGraph

import (
	"log"
	"math"
	"sync"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

func SmithWaterman(alpha []dna.Base, beta []dna.Base, config *GraphSettings, pool *sync.Pool) (int64, []cigar.Cigar, int, int, int, int) {
	var rows, columns int = len(alpha) + 1, len(beta) + 1
	var i, j, minI, maxI, minJ, maxJ int
	var currMax int64 = math.MinInt64

	dp := pool.Get().(*Matrix)
	defer pool.Put(dp)
	dp.Reset(rows, columns)

	for i = 1; i < rows; i++ {
		for j = 1; j < columns; j++ {
			dp.matrix[i][j], dp.trace[i][j] = cigar.TripleMaxTraceExtended(dp.matrix[i-1][j-1], dp.matrix[i-1][j-1]+config.ScoreMatrix[alpha[i-1]][beta[j-1]], dp.matrix[i][j-1]+config.GapPenalty, dp.matrix[i-1][j]+config.GapPenalty)
			if dp.matrix[i][j] > currMax {
				currMax = dp.matrix[i][j]
				maxI = i
				maxJ = j
			}
			if dp.matrix[i][j] < 0 {
				dp.matrix[i][j] = 0
			}
		}
	}
	route := make([]cigar.Cigar, 1)
	route[0] = cigar.Cigar{RunLength: 0, Op: dp.trace[maxI][maxJ]}

	for i, j = maxI, maxJ; dp.matrix[i][j] > 0; {
		if route[len(route)-1].Op != dp.trace[i][j] {
			route = append(route, cigar.Cigar{RunLength: 1, Op: dp.trace[i][j]})
		} else {
			route[len(route)-1].RunLength++
		}
		switch dp.trace[i][j] {
		case cigar.Equal, cigar.Mismatch:
			i, j = i-1, j-1
		case cigar.Insertion:
			j--
		case cigar.Deletion:
			i--
		default:
			log.Fatalf("Error: unexpected traceback: %c", dp.trace[i][j])
		}
		minI, minJ = i, j
	}
	cigar.ReverseCigar(route)
	return dp.matrix[maxI][maxJ], route, minI, maxI, minJ, maxJ
}

func LeftLocal(alpha []dna.Base, beta []dna.Base, config *GraphSettings, pool *sync.Pool) (int64, []cigar.Cigar, int, int, int, int) {
	var i, j int
	rows, columns := len(alpha), len(beta)

	dp := pool.Get().(*Matrix)
	defer pool.Put(dp)
	dp.Reset(rows+1, columns+1)

	for i = 1; i <= rows; i++ {
		for j = 1; j <= columns; j++ {
			dp.matrix[i][j], dp.trace[i][j] = cigar.TripleMaxTraceExtended(dp.matrix[i-1][j-1], dp.matrix[i-1][j-1]+config.ScoreMatrix[alpha[i-1]][beta[j-1]], dp.matrix[i][j-1]+config.GapPenalty, dp.matrix[i-1][j]+config.GapPenalty)
			if dp.matrix[i][j] < 0 {
				dp.matrix[i][j] = 0
			}
		}
	}
	var minI, minJ = rows, columns
	route := make([]cigar.Cigar, 1)
	route[0] = cigar.Cigar{RunLength: 0, Op: dp.trace[minI][minJ]}
	for i, j = minI, minJ; dp.matrix[i][j] > 0; {
		if route[len(route)-1].Op != dp.trace[i][j] {
			route = append(route, cigar.Cigar{RunLength: 1, Op: dp.trace[i][j]})
		} else {
			route[len(route)-1].RunLength++
		}
		switch dp.trace[i][j] {
		case cigar.Equal, cigar.Mismatch:
			i, j = i-1, j-1
		case cigar.Insertion:
			j--
		case cigar.Deletion:
			i--
		default:
			log.Fatalf("Error: unexpected traceback %c\n", dp.trace[i][j])
		}
		minI = i
		minJ = j
	}
	cigar.ReverseCigar(route)
	return dp.matrix[len(alpha)][len(beta)], route, minI, rows, minJ, columns
}

func RightLocal(alpha []dna.Base, beta []dna.Base, config *GraphSettings, pool *sync.Pool) (int64, []cigar.Cigar, int, int, int, int) {
	rows, columns := len(alpha)+1, len(beta)+1
	dp := pool.Get().(*Matrix)
	defer pool.Put(dp)
	dp.Reset(rows, columns)

	var currMax int64 = math.MinInt64
	var i, j, maxI, maxJ int

	for i = 0; i < rows; i++ {
		for j = 0; j < columns; j++ {
			if i == 0 && j == 0 {
				dp.matrix[i][j] = 0
			} else if i == 0 {
				dp.matrix[i][j] = dp.matrix[i][j-1] + config.GapPenalty
				dp.trace[i][j] = cigar.Insertion
			} else if j == 0 {
				dp.matrix[i][j] = dp.matrix[i-1][j] + config.GapPenalty
				dp.trace[i][j] = cigar.Deletion
			} else {
				dp.matrix[i][j], dp.trace[i][j] = cigar.TripleMaxTraceExtended(dp.matrix[i-1][j-1], dp.matrix[i-1][j-1]+config.ScoreMatrix[alpha[i-1]][beta[j-1]], dp.matrix[i][j-1]+config.GapPenalty, dp.matrix[i-1][j]+config.GapPenalty)
			}
			if dp.matrix[i][j] > currMax {
				currMax = dp.matrix[i][j]
				maxI = i
				maxJ = j
			}
		}
	}
	route := make([]cigar.Cigar, 1)
	route[0] = cigar.Cigar{RunLength: 0, Op: dp.trace[maxI][maxJ]}
	for i, j = maxI, maxJ; i > 0 || j > 0; {
		if route[len(route)-1].Op != dp.trace[i][j] {
			route = append(route, cigar.Cigar{RunLength: 1, Op: dp.trace[i][j]})
		} else {
			route[len(route)-1].RunLength++
		}
		switch dp.trace[i][j] {
		case cigar.Equal, cigar.Mismatch:
			i, j = i-1, j-1
		case cigar.Insertion:
			j--
		case cigar.Deletion:
			i--
		default:
			log.Fatalf("Error: unexpected traceback with %c\n", dp.trace[i][j])
		}
	}
	cigar.ReverseCigar(route)
	return dp.matrix[maxI][maxJ], route, 0, maxI, 0, maxJ
}
