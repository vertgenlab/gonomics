package align

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

func ConstGap_lowMem(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64) (int64, []Cigar) {
	//make matrices for m (score, broken down into mRowCurrent,mRowPrevious and mColumn to save memory) and trace (traceback with directions)
	mRowCurrent := make([]int64,len(beta)+1)
	mRowPrevious := make([]int64,len(beta)+1)
	var mColumn int = len(alpha)+1
	trace := make([][]ColType, len(alpha)+1)
	for idx := 0; idx < mColumn; idx++ {
		trace[idx] = make([]ColType, len(beta)+1)
	}

	//fill out matrices for m and trace
	var i, j, routeIdx int
	for i = 0; i < mColumn; i++ {
		for j = range mRowCurrent {
			if i == 0 && j == 0 {
				mRowCurrent[j] = 0
			} else if i == 0 {
				mRowCurrent[j] = mRowCurrent[j-1] + gapPen
				trace[i][j] = 1
			} else if j == 0 {
				mRowCurrent[j] = mRowPrevious[j] + gapPen
				trace[i][j] = 2
			} else {
				mRowCurrent[j], trace[i][j] = tripleMaxTrace(mRowPrevious[j-1]+scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1]+gapPen, mRowPrevious[j]+gapPen)
			}
		}
		if i < mColumn-1 {
			mRowPrevious,mRowCurrent=mRowCurrent,mRowPrevious  //reuse mRow variables to save memory, swap slices rather than copy(mRowPrevious,mRowCurrent) for easy update, but only up until the second to last row
		}
	}

	//write cigar
	route := make([]Cigar, 1)
	for i, j, routeIdx = mColumn-1, len(mRowCurrent)-1, 0; i > 0 || j > 0; {
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
	return mRowCurrent[len(mRowCurrent)-1], route
}
