package align

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"math"
)

func ConstGap(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64) (int64, []Cigar) {

	var checkersize int
	checkersize = 10000 //TODO: make checkersize an cmd option; checkersize_i and checkersize_j

	//Step 1: find highest score, as well as get the position (i and j) of the highest score, and materials needed to fill traceback and write cigar in checkerboards
	score_highest, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j := HighestScore(alpha, beta, scores, gapPen, checkersize)

	//Make variables needed for checkerboards
	//k1: the i-index (row-index) of the current checkerboard
	//k2: the j-index (column-index) of the current checkerboard
	//i_inChecker_max: the max i-index (row-index) in the current checkerboard during Step 2 FillTraceback, either checkersize-1 or a smaller number, e.g. if the alpha/beta sequence to be aligned is not divisible into perfect checkerboards, and the remainder checkerboard is not completely filled
	//j_inChecker_max: the max j-index (column-index) in the current checkerboard during Step 2 FillTraeback
	//i_inChecker_min: the min i-index (row-index) in the current checkerboard during Step 3 WriteCigar, either 0 or another number, e.g. if cigar route leaves the checkerboard at a position that is not i_inChecker==0, j_inChecker==0
	//j_inChecker_min: the min j-index (column-index) in the current checkerboard during Step 3 WriteCigar
	var k1, k2, i_inChecker_max, j_inChecker_max, i_inChecker_min, j_inChecker_min int
	i_inChecker_min = -2                    //initialize i_inChecker_min != 0, so that the first ever Step 3 WriteCigar will not interfere with i_inChecker_max
	j_inChecker_min = -2                    //ditto for j
	trace := make([][]ColType, checkersize) //trace matrix size is checkersize*checkersize
	for idx := 0; idx < len(trace); idx++ {
		trace[idx] = make([]ColType, checkersize)
	}
	route := make([]Cigar, 1) //initialie cigar route and routeIdx
	var routeIdx_current int = 0

	for k1, k2 = int((score_highest_i-1)/checkersize), int((score_highest_j-1)/checkersize); k1 >= 0 && k2 >= 0; { //use a function of score_highest_i, score_highest_j, and checkersize to initialize the right k1 and k2, go to the correct checkerboard to start traceback

		//Step 2: for a checkerboard, fill traceback, as well as find the max i and j when the checkerboard traceback ends
		trace, i_inChecker_max, j_inChecker_max = FillTraceback(alpha, beta, scores, gapPen, checkersize, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j, k1, k2, i_inChecker_min, j_inChecker_min)

		//Step 3: for a checkerboard, use traceback to write cigar (update route and routeIdx), as well as find the min i and j when the checkerboard cigar ends
		route, routeIdx_current, i_inChecker_min, j_inChecker_min = WriteCigar(trace, i_inChecker_max, j_inChecker_max, route, routeIdx_current, i_inChecker_min, j_inChecker_min)

		//Use Step 3's i_inChecker_min and j_inChecker_min to find the next checkerboard to go to and where to start in that checkerboard (update k1, k2)
		if i_inChecker_min < 0 && j_inChecker_min < 0 {
			k1, k2 = k1-1, k2-1 //go to the next checkerboard with lower i, lower j
		} else if i_inChecker_min < 0 { //but j_inChecker != 0
			k1 -= 1 //go to the next checkerboard with lower i, same j
		} else if j_inChecker_min < 0 { //but i_inChecker != 0
			k2 -= 1 //go to the next checkerboard with lower j, same i
		}

	}

	//Step 4: write the last cigar entry
	//This step is necessary because the row i=0 and the column j=0 are always stored in trace_prep, and never filled in in any checkerboard, but the cigar ends by the route going into 1 box/entry in either the row i=0 or the column j=0
	if i_inChecker_min != -1 && j_inChecker_min == -1 { //indicating that Step 3 ended when k1 and k2 reached the smallest combination, and reached j=0, so the last cigar is a "2", e.g. if last cigar entry is the i=1, j=0 square
		route, routeIdx_current = LastCigar(route, routeIdx_current, 2)
	} else if i_inChecker_min == -1 && j_inChecker_min != -1 { //indicating that Step 3 ended when k1 and k2 reached the smallest combination, and reached i=0, so the last cigar is a "1", e.g. if last cigar entry is the i=0, j=1 square
		route, routeIdx_current = LastCigar(route, routeIdx_current, 1)
	} //no more "else" because the only situation left is if Step 3 ended when k1 and k2 reached the smallest combination, and reached both i=0 and j=0, aka the i=0 j=0 square, and there is no sequence there, so no cigar

	//Final processing (reverse route) and return outputs
	reverseCigar(route)
	return score_highest, route
}

//Step 1
func HighestScore(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, checkersize int) (int64, int, int, [][]int64, [][]int64) {

	mRowCurrent := make([]int64, len(beta)+1)
	mRowPrevious := make([]int64, len(beta)+1)
	var mColumn int = len(alpha) + 1
	trace_prep_i := make([][]int64, int(len(alpha)/checkersize)+1) //trace_prep_i saves all rows (i) that are needed to initialize checkerboard tracebacks
	trace_prep_j := make([][]int64, int(len(beta)/checkersize)+1)  //trace_prep_j saves all columns (j) that are needed to initialize checkerboard tracebacks
	for idx := 0; idx < len(trace_prep_i); idx++ {
		trace_prep_i[idx] = make([]int64, len(beta)+1)
	}
	for idx := 0; idx < len(trace_prep_j); idx++ {
		trace_prep_j[idx] = make([]int64, len(alpha)+1)
	}
	var i, j int

	for i = 0; i < mColumn; i++ {

		for j = range mRowCurrent {
			if i == 0 && j == 0 {
				mRowCurrent[j] = 0
				trace_prep_j[j/checkersize][i] = mRowCurrent[j] //it is implied here that j%checkersize==0. It must be saved in trace_prep_j
			} else if i == 0 {
				mRowCurrent[j] = mRowCurrent[j-1] + gapPen
				if j%checkersize == 0 {
					trace_prep_j[j/checkersize][i] = mRowCurrent[j]
				}
			} else if j == 0 {
				mRowCurrent[j] = mRowPrevious[j] + gapPen
				trace_prep_j[j/checkersize][i] = mRowCurrent[j] //it is implied that j%checkersize==0. It must be saved in trace_prep_j
			} else {
				mRowCurrent[j], _ = tripleMaxTrace(mRowPrevious[j-1]+scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1]+gapPen, mRowPrevious[j]+gapPen)
				if j%checkersize == 0 {
					trace_prep_j[j/checkersize][i] = mRowCurrent[j]
				}
			}
		}

		if i%checkersize == 0 && i < mColumn-1 {
			copy(trace_prep_i[i/checkersize], mRowCurrent) //copy instead of assign values to avoid pointers automatically updating/synchronizing the 2 matrices
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
		} else if i < mColumn-1 {
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious //reuse mRow variables to save memory, but only up until the second to last row
		}
	}

	score_highest := mRowCurrent[len(mRowCurrent)-1]
	score_highest_i := mColumn - 1
	score_highest_j := len(mRowCurrent) - 1
	return score_highest, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j
}

//Step 2
func FillTraceback(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, checkersize int, score_highest_i int, score_highest_j int, trace_prep_i [][]int64, trace_prep_j [][]int64, k1 int, k2 int, i_inChecker_min_Previous int, j_inChecker_min_Previous int) ([][]ColType, int, int) {

	mRowCurrent := make([]int64, len(beta)+1)
	mRowPrevious := make([]int64, len(beta)+1)
	mRowPrevious = trace_prep_i[k1]
	//i_inChecker and j_inChecker are internal (within checkerboard) versions of i and j, e.g. for checkersize=3, i_inChecker can take on the values 0, 1, 2
	//i_max and j_max are external (non-inChecker) versions of i_inChecker_max and j_inChecker_max
	var i, j, i_inChecker, j_inChecker, i_inChecker_max, j_inChecker_max, i_max, j_max int
	if i_inChecker_min_Previous >= 0 {
		//use the i_inChecker_min_Previous from the previous k1/k2 loop iteration's Step 3 WriteCigar to find the i_max at which to stop filling the trace matrix
		i_max = checkersize*k1 + 1 + i_inChecker_min_Previous
	} else {
		//otherwise, since no restriction from previous k1/k2 loop Step 3 WriteCigar, find the minimum value between "where the next checkerboard starts" and "the i of the highest score", and that is where to stop filing the trace matrix
		i_max = int(math.Min(float64(checkersize*(k1+1)), float64(score_highest_i)))
	}
	if j_inChecker_min_Previous >= 0 {
		j_max = checkersize*k2 + 1 + j_inChecker_min_Previous
	} else {
		j_max = int(math.Min(float64(checkersize*(k2+1)), float64(score_highest_j)))
	}
	i_inChecker_max = (i_max - 1) % checkersize
	j_inChecker_max = (j_max - 1) % checkersize
	trace := make([][]ColType, checkersize)
	for idx := 0; idx < len(trace); idx++ {
		trace[idx] = make([]ColType, checkersize)
	}

	for i = checkersize*k1 + 1; i <= i_max; i++ {
		i_inChecker = (i - 1) % checkersize
		mRowCurrent[checkersize*k2] = trace_prep_j[k2][checkersize*k1+1+i_inChecker]

		for j = checkersize*k2 + 1; j <= j_max; j++ {
			j_inChecker = (j - 1) % checkersize
			mRowCurrent[j], trace[i_inChecker][j_inChecker] = tripleMaxTrace(mRowPrevious[j-1]+scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1]+gapPen, mRowPrevious[j]+gapPen) //it is ok even if mRowCurrent isn't completely filled, aka only part of mRowCurrent is needed for the checkerboard traceback and writing cigar
		}

		if i <= checkersize*(k1+1)-1 && i <= score_highest_i-1 {
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
		}
	}

	return trace, i_inChecker_max, j_inChecker_max
}

//Step 3
func WriteCigar(trace [][]ColType, i_inChecker_max_Previous int, j_inChecker_max_Previous int, route []Cigar, routeIdx_current int, i_inChecker_min_Previous int, j_inChecker_min_Previous int) ([]Cigar, int, int, int) {

	var i_inChecker, j_inChecker, routeIdx, i_inChecker_min, j_inChecker_min, i_inChecker_max, j_inChecker_max int
	route_updated := route

	//use the i_inChecker_min_Previous from the previous k1/k2 loop iteration's Step 3 WriteCigar to find the i_inChecker_max at which to start writing cigar
	if i_inChecker_min_Previous >= 0 {
		i_inChecker_max = i_inChecker_min_Previous
	} else {
		i_inChecker_max = i_inChecker_max_Previous
	}
	if j_inChecker_min_Previous >= 0 {
		j_inChecker_max = j_inChecker_min_Previous
	} else {
		j_inChecker_max = j_inChecker_max_Previous
	}

	for i_inChecker, j_inChecker, routeIdx = i_inChecker_max, j_inChecker_max, routeIdx_current; i_inChecker >= 0 && j_inChecker >= 0; {

		//write cigar segment
		if route_updated[routeIdx].RunLength == 0 {
			route_updated[routeIdx].RunLength = 1
			route_updated[routeIdx].Op = trace[i_inChecker][j_inChecker]
		} else if route_updated[routeIdx].Op == trace[i_inChecker][j_inChecker] {
			route_updated[routeIdx].RunLength += 1
		} else {
			route_updated = append(route_updated, Cigar{RunLength: 1, Op: trace[i_inChecker][j_inChecker]})
			routeIdx++
		}

		//update i_inChecker and j_inChecker
		switch trace[i_inChecker][j_inChecker] {
		case 0:
			i_inChecker, j_inChecker = i_inChecker-1, j_inChecker-1
		case 1:
			j_inChecker -= 1
		case 2:
			i_inChecker -= 1
		default:
			log.Fatalf("Error: unexpected traceback")
		}

		i_inChecker_min = i_inChecker
		j_inChecker_min = j_inChecker

	}

	return route_updated, routeIdx, i_inChecker_min, j_inChecker_min //return the routeIdx that was updated in the loop, not the function input "routeIdx_current"
}

//Step 4
func LastCigar(route []Cigar, routeIdx_current int, Op_end ColType) ([]Cigar, int) {

	if route[routeIdx_current].Op == Op_end {
		route[routeIdx_current].RunLength += 1
	} else {
		route = append(route, Cigar{RunLength: 1, Op: Op_end})
		routeIdx_current++
	}

	return route, routeIdx_current
}
