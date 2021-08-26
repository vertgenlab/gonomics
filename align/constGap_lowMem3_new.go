package align

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"fmt" //TODO: remove after debugging
)

func ConstGap(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64) (int64, []Cigar) {
	var checkersize int
	checkersize = 3 //TODO: adjust hardcode after debugging. Default to 10000. Also maybe consider different checkersize_i and checkersize_j

	//Step 1
	score_highest, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j := HighestScore(alpha, beta, scores, gapPen, checkersize)

	//Steps 2 an 3 are in a loop, start Step 2 with highest checkerboard coordinates, end at Step 3 when cigar reaches i==0 or j==0
	var k1, k2 int
	var i_inChecker_max, j_inChecker_max int
	var i_inChecker_min, j_inChecker_min int
	trace := make([][]ColType, checkersize)
	for idx := 0; idx < len(trace); idx++ {
		trace[idx] = make([]ColType, checkersize)
	}
	route := make([]Cigar, 1)
	var routeIdx_max int = 0

	for k1, k2 = len(trace_prep_i)-1, len(trace_prep_j)-1; k1 >= 0 || k2 >= 0; { //while loop. TODO: is the end condition right? Is this the end of traceback? Worked for test example1. Also, replace initial k1 and k2 with int(score_highest_i/checkersize)
		//Try to fix test2 by skipping when last lines are multiples of checkersize
		//if checkersize*k1+1 != score_highest_i { //TODO: should delete when make sure it's not needed. Pretty sure it's not needed
		//Step 2
		trace, i_inChecker_max, j_inChecker_max = FillTraceback(alpha, beta, scores, gapPen, checkersize, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j, k1, k2)

		//Try updating k1, k2 based on i_inChecker_max too, not just i_inChecker_min like below
		if j_inChecker_max == 0 { //it's impossible that both i and j inChecker max are 0, because then there won't be a checkerboard. j_inChecker_max==0 is the sitaution in test2, need to directly update k because can't get useful cigar
			// go to the next checkerboard with lower j, same i
			k2 -= 1
		//} else if i_inChecker_max == 0 { //but j_inChecker_max != 0, but this doesn't seem to matter based on how I wrote mRowCurrent, like this can still be put into WriteCigar and get useful cigar
			// go to the next checkerboard with lower i, same j
		//	k1 -= 1
		} else {

			//Step 3
			route, routeIdx_max, i_inChecker_min, j_inChecker_min = WriteCigar(trace, i_inChecker_max, j_inChecker_max, route, routeIdx_max)
			//}
			//Update k1,k2
			if i_inChecker_min == 0 && j_inChecker_min == 0 { //when cigar is done, find out how while loop ended, and update k1, k2
				// go to the next checkerboard with lower i, lower j
				k1, k2 = k1-1, k2-1
			} else if i_inChecker_min == 0 { //but j_inChecker != 0
				// go to the next checkerboard with lower i, same j
				k1 -= 1
			} else if j_inChecker_min == 0 { //but i_inChecker != 0
				// go to the next checkerboard with lower j, same i
				k2 -= 1
			}
			//TODO: add another feature, instead of always filling up an entire checkerboard, fill up to the highest/rightest box indicated by the last checkerboard's i_inChecker_min and j_inChecker_min
		}
	}

	//After Step 3 ends
	reverseCigar(route)
	fmt.Printf("ConstGap final route: %v\n", route)
	return score_highest, route
}

//Step 1
func HighestScore(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, checkersize int) (int64, int, int, [][]int64, [][]int64) {
	mRowCurrent := make([]int64, len(beta)+1)
	mRowPrevious := make([]int64, len(beta)+1)
	var mColumn int = len(alpha) + 1
	trace_prep_i := make([][]int64,int(len(alpha)/checkersize)+1)
	trace_prep_j := make([][]int64,int(len(beta)/checkersize)+1)
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
				trace_prep_j[j/checkersize][i] = mRowCurrent[j] //implied that j%checkersize==0
			} else if i == 0 {
				mRowCurrent[j] = mRowCurrent[j-1] + gapPen
				if j % checkersize == 0 {
					trace_prep_j[j/checkersize][i] = mRowCurrent[j]
				}
			} else if j == 0 {
				mRowCurrent[j] = mRowPrevious[j] + gapPen
				trace_prep_j[j/checkersize][i] = mRowCurrent[j] //implied that j%checkersize==0
			} else {
				mRowCurrent[j], _ = tripleMaxTrace(mRowPrevious[j-1]+scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1]+gapPen, mRowPrevious[j]+gapPen)
				if j % checkersize == 0 {
					trace_prep_j[j/checkersize][i] = mRowCurrent[j]
				}
			}
		}
		fmt.Printf("Before fill: i, i_remainder_checkersize, mRowCurrent, trace_prep_i : %d, %d, %v, %v\n", i, i%checkersize, mRowCurrent, trace_prep_i) //TODO: why does trace_prep_i change before fill? are ifs executed in parallel? is "=" not one-time? Maybe pointer memory operation in the background.
		if i % checkersize == 0 && i < mColumn-1 {
			copy(trace_prep_i[i/checkersize],mRowCurrent) //by coping value, instead of assigning "trace_prep_i[i/checkersize] = mRowCurrent", solves trace_prep_i problem
			fmt.Printf("After fill: i, i_remainder_checkersize, mRowCurrent, trace_prep_i : %d, %d, %v, %v\n", i, i%checkersize, mRowCurrent, trace_prep_i) //TODO
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
		} else if i < mColumn-1 {
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious //reuse mRow variables to save memory, but only up until the second to last row
		}
	}
	score_highest := mRowCurrent[len(mRowCurrent)-1]
	score_highest_i := mColumn-1
	score_highest_j := len(mRowCurrent)-1
	fmt.Printf("score_highest: %d\n", score_highest) //TODO: remove after debugging. PASSED
	fmt.Printf("trace_prep_i: %v\n", trace_prep_i) //TODO: remove after debugging. PASSED
	fmt.Printf("trace_prep_j: %v\n", trace_prep_j) //TODO: remove after debugging. PASSED

	return score_highest, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j
}

//Step 2
func FillTraceback(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64, checkersize int, score_highest_i int, score_highest_j int, trace_prep_i [][]int64, trace_prep_j [][]int64, k1 int, k2 int) ([][]ColType, int, int) {
	mRowCurrent := make([]int64, len(beta)+1)
	mRowPrevious := make([]int64, len(beta)+1)
	mRowPrevious = trace_prep_i[k1]
	var i, j int
	var i_inChecker, j_inChecker int
	var i_inChecker_max, j_inChecker_max int //to keep track of i and j max in a checkerboard
	trace := make([][]ColType, checkersize)
	for idx := 0; idx < len(trace); idx++ {
		trace[idx] = make([]ColType, checkersize)
	}

	for i = checkersize*k1+1; i <= checkersize*(k1+1) && i <= score_highest_i; i++ {
		i_inChecker = (i-1) % checkersize
		i_inChecker_max = i_inChecker //update i_inChecker_max. TODO: can I, and is it better to get this data outside of loop?
		mRowCurrent[checkersize*k2] = trace_prep_j[k2][checkersize*k1+1+i_inChecker]
		fmt.Printf("initialize checkerboard: k1, k2, i, trace_prep_i[k1], mRowPrevious, mRowCurrent: %d, %d, %d, %v, %v, %v\n", k1, k2, i, trace_prep_i[k1], mRowPrevious, mRowCurrent) //TODO

		for j = checkersize*k2+1; j <= checkersize*(k2+1) && j <= score_highest_j; j++ {
			j_inChecker = (j-1) % checkersize
			j_inChecker_max = j_inChecker
			mRowCurrent[j], trace[i_inChecker][j_inChecker] = tripleMaxTrace(mRowPrevious[j-1]+scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1]+gapPen, mRowPrevious[j]+gapPen) //it is ok if mRowCurrent isn't completely filled, aka only part of mRowCurrent is needed for the checkerboard
			fmt.Printf("filled checkerboard: k1, k2, i, j, i_inChecker, j_inChecker, alpha[i-1], beta[j-1], mRowCurrent[j], trace[i_inChecker][j_inChecker]: %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", k1, k2, i, j, i_inChecker, j_inChecker, alpha[i-1], beta[j-1], mRowCurrent[j], trace[i_inChecker][j_inChecker]) //TODO: remove after debugging
			//fmt.Printf("mRowPrevious[j-1], scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1], gapPen, mRowPrevious[j]: %d, %d, %d, %d, %d\n", mRowPrevious[j-1], scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1], gapPen, mRowPrevious[j])
		}
		if i <= checkersize*(k1+1)-1 && i <= score_highest_i-1 {
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
			fmt.Printf("switched\n")
		}
	}
	return trace, i_inChecker_max, j_inChecker_max
}

//Step 3: given trace, i_inChecker_max, j_inChecker_max (can be smaller than checkersize even in a central checkerboard), route, routeIdx, return cigar (route, routeIdx), next checkerboard coordinates (final i_inChecker_min, j_inChecker_min)
func WriteCigar(trace [][]ColType, i_inChecker_max int, j_inChecker_max int, route []Cigar, routeIdx_max int) ([]Cigar, int, int, int) {
	var i_inChecker, j_inChecker, routeIdx int
	var i_inChecker_min, j_inChecker_min int
	route_updated := route

	fmt.Printf("entered WriteCigar: i_inChecker_max, j_inChecker_max, routeIdx_max: %d, %d, %d\n", i_inChecker_max, j_inChecker_max, routeIdx_max)
	for i_inChecker, j_inChecker, routeIdx = i_inChecker_max, j_inChecker_max, routeIdx_max; i_inChecker >= 0 && j_inChecker >= 0; { //Did this change lead to infinite loop? No

		i_inChecker_min = i_inChecker //TODO: can I, and is it better to get this data outside of loop?/inside of loop?
		j_inChecker_min = j_inChecker

		if route_updated[routeIdx].RunLength == 0 {
			route_updated[routeIdx].RunLength = 1
			route_updated[routeIdx].Op = trace[i_inChecker][j_inChecker]
		} else if route_updated[routeIdx].Op == trace[i_inChecker][j_inChecker] {
			route_updated[routeIdx].RunLength += 1
		} else {
			route_updated = append(route_updated, Cigar{RunLength: 1, Op: trace[i_inChecker][j_inChecker]})
			routeIdx++
			fmt.Printf("WriteCigar while loop: i_inChecker, j_inChecker, routeIdx, route_updated: %d, %d, %d, %v\n", i_inChecker, j_inChecker, routeIdx, route_updated)
		}

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
	}

	fmt.Printf("about to leave WriteCigar: i_inChecker_min, j_inChecker_min, routeIdx, route_updated: %d, %d, %d, %v\n", i_inChecker_min, j_inChecker_min, routeIdx, route_updated)

	return route_updated, routeIdx, i_inChecker_min, j_inChecker_min //return the routeIdx that was updated in the loop, not the input routeIdx_max
}
