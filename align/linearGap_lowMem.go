package align

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

func ConstGap(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64) (int64, []Cigar) {
	//make matrices for m (score) and trace (traceback with directions)
	//m := make([][]int64, len(alpha)+1)
	mRowCurrent := make([]int64,len(beta)+1)
	mRowPrevious := make([]int64,len(beta)+1)
	//var mColumn int = len(alpha)+1
	trace := make([][]ColType, len(alpha)+1)
	//for idx := range m {
	for idx := range trace {
		//m[idx] = make([]int64, len(beta)+1)
		trace[idx] = make([]ColType, len(beta)+1)
	}

	//fill out matrices for m and trace
	var i, j, routeIdx int
	for i = 0; i < len(alpha)+1; i++ {
		for j = 0; j < len(beta)+1; j++ {
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
		if i < len(alpha) {
			mRowPrevious,mRowCurrent=mRowCurrent,mRowPrevious
		} //swap slices rather than copy(mRowPrevious,mRowCurrent) for easy update, but only up until the second to last row
	}
	//TODO: cleanup for loop
	//for j = 0; j < mColumn; j++ {
	//for j = range m[0] {
	//raven's note: by getting rid of 1 for loop, I now have 1 fewer layer of {}
			//if i == 0 && j == 0 {
	//		if mRowCurrent[j] == 0 && j == 0 {
	//			mRowCurrent[j] = 0
			//} else if i == 0 {
	//		} else if mRowCurrent[j] == 0 { //raven's note: but j!=0
				//m[i][j] = m[i][j-1] + gapPen
	//			mRowCurrent[j] = mRowCurrent[j-1] + gapPen
	//			trace[i][j] = 1
	//		} else if j == 0 { //raven's note: but mRowCurrent[j]!=0
				//m[i][j] = m[i-1][j] + gapPen
	//			mRowCurrent[j] = mRowPrevious[j] + gapPen
	//			trace[i][j] = 2
	//		} else { //raven's note: mRowCurrent[j]!=0 && j!=0
				//m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
	//			mRowCurrent[j], trace[i][j] = tripleMaxTrace(mRowPrevious[j-1]+scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1]+gapPen, mRowPrevious[j]+gapPen)
	//		}
	//		copy(mRowPrevious,mRowCurrent) //raven's note: refresh previous so it takes on the value of current at the end of each row. Faster to swap them "a,b=b,a"
	//}

	//write cigar
	route := make([]Cigar, 1)
	//for i, j, routeIdx = len(trace)-1, len(trace[0])-1, 0; i > 0 || j > 0; {
	//raven's note for above code: len(trace)=len(beta)+1
	for i, j, routeIdx = len(alpha), len(beta), 0; i > 0 || j > 0; {
	//for i, j, routeIdx = len(beta), len(alpha), 0; i > 0 || j > 0; { //raven's note: is it for i=len(alpha) or i=len(beta)?
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
			//i -= 1
		case 2:
			i -= 1
			//j -= 1
		default:
			log.Fatalf("Error: unexpected traceback")
		}
	}
	reverseCigar(route)
	//return m[len(m)-1][len(m[0])-1], route //m[len(m)-1][len(m[0]-1)]==m[len(alpha)][len(beta)]
	return mRowCurrent[len(mRowCurrent)-1], route
}
