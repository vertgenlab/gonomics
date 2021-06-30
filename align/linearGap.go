package align

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

func ConstGap_testing(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64) (int64, []Cigar) {
	m := make([][]int64, len(alpha)+1)
	trace := make([][]ColType, len(alpha)+1)
	for idx := range m {
		m[idx] = make([]int64, len(beta)+1)
		trace[idx] = make([]ColType, len(beta)+1)
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
				m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			}
		}
	}

	route := make([]Cigar, 1)
	for i, j, routeIdx = len(trace)-1, len(trace[0])-1, 0; i > 0 || j > 0; {
		if route[routeIdx].RunLength == 0 {
			route[routeIdx].RunLength = 1
			route[routeIdx].Op = trace[i][j]
		} else if route[routeIdx].Op == trace[i][j] {
			route[routeIdx].RunLength += 1
		} else {
			route = append(route, Cigar{RunLength: 1, Op: trace[i][j]})
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
