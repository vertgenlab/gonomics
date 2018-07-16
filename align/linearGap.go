package align

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
)

func ConstGap(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64) (int64, []cigar) {
	m := make([][]int64, len(alpha)+1)
	trace := make([][]colType, len(alpha)+1)
	for idx := range m {
		m[idx] = make([]int64, len(beta)+1)
		trace[idx] = make([]colType, len(beta)+1)
	}

	var i, j, routeIdx int
	for i, _ = range m {
		for j, _ = range m[0] {
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

	route := make([]cigar, 1)
	for i, j, routeIdx = len(trace)-1, len(trace[0])-1, 0; i > 0 || j > 0; {
		if route[routeIdx].runLength == 0 {
			route[routeIdx].runLength = 1
			route[routeIdx].op = trace[i][j]
		} else if route[routeIdx].op == trace[i][j] {
			route[routeIdx].runLength += 1
		} else {
			route = append(route, cigar{runLength: 1, op: trace[i][j]})
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
			common.Exit("Error: unexpected traceback")
		}
	}
	reverseCigar(route)
	return m[len(m)-1][len(m[0])-1], route
}
