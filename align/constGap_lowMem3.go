package align

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"fmt" //TODO: remove after debugging
)

//e.g. constGap_alignTests_8. len(alpha)=10,len(beta)=14
func ConstGap_testing3(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapPen int64) (int64, []Cigar) {
	//make matrices for m (score, broken down into mRowCurrent,mRowPrevious and mColumn to save memory)
	//and trace (traceback with directions. Make prep matrix for every checkersize, and trace matrix representing checkersize * len(beta) box)
	var checkersize int
	checkersize = 3 //TODO: adjust hardcode after debugging. Default to 10000. Also maybe consider different checkersize_i and checkersize_j
	mRowCurrent := make([]int64, len(beta)+1) //e.g. len(beta)+1=15
	mRowPrevious := make([]int64, len(beta)+1)
	var mColumn int = len(alpha) + 1 //len(alpha)+1=11
	//divide len(alpha)/checkersize and take int answer, no remainder,+1, that is how many rows are needed in the trace_prep_i matrix
	trace_prep_i := make([][]int64,int(len(alpha)/checkersize)+1) //e.g. should save rows i=0,3,6,9, 4 rows, 10/3+1=4 (+1 for i=0)
	//similarly, need trace_prep_j matrix, but unlike trace_prep_i, for trace_prep_j the outer layer is column, and the inner layer is row
	trace_prep_j := make([][]int64,int(len(beta)/checkersize)+1) //e.g. should save columns j=0,3,6,9,12, 5 columns, 14/3+1=5 (+1 for j=0)
	for idx := 0; idx < len(trace_prep_i); idx++ {
		//save the entire row in trace_prep_i, parse out those positiosn needed for checkersize later
		trace_prep_i[idx] = make([]int64, len(beta)+1) //e.g. trace_prep_j is 5 columns * 11 rows
	}
	for idx := 0; idx < len(trace_prep_j); idx++ {
		trace_prep_j[idx] = make([]int64, len(alpha)+1)
	}
	trace := make([][]ColType, checkersize) //trace matrix is of size checkersize
	for idx := 0; idx < len(trace); idx++ {
		trace[idx] = make([]ColType, checkersize) //e.g. trace is 3 rows * 3 columns
	}

	//get highest score. TODO: make this its own function return highest score (including position ij, not necessarily top right corner if not forcing ends to align) and trace_prep_i, trace_prep_j
	//Make another function that can do 1 checkerboard square and determine which square to go to next (including position, which j in mRowCurrent)
	var i, j int
	for i = 0; i < mColumn; i++ {
		for j = range mRowCurrent {
			if i == 0 && j == 0 {
				mRowCurrent[j] = 0
				trace_prep_j[j/checkersize][i] = mRowCurrent[j] //for sure j%checkersize ==0, aka trace_prep_j[0][0]=0
			} else if i == 0 {
				mRowCurrent[j] = mRowCurrent[j-1] + gapPen
				if j % checkersize == 0 {
					trace_prep_j[j/checkersize][i] = mRowCurrent[j] //e.g. j=3,6,9,12, i=0, get those scores into the trace_prep_j[1,2,3,4][0]
				}
			} else if j == 0 {
				mRowCurrent[j] = mRowPrevious[j] + gapPen
				trace_prep_j[j/checkersize][i] = mRowCurrent[j] //for sure j%checkersize ==0, e.g. get j=0, i=all except 0 into trace_prep_j[0][all except 0]
			} else {
				mRowCurrent[j], _ = tripleMaxTrace(mRowPrevious[j-1]+scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1]+gapPen, mRowPrevious[j]+gapPen)
				if j % checkersize == 0 {
					trace_prep_j[j/checkersize][i] = mRowCurrent[j] //e.g. j=3,6,9,12, i=all except 0, get those scores into the trace_prep_j[1,2,3,4][all except 0]. By this point trace_prep_j is all filled
				}
			}
		}
		fmt.Printf("Before fill: i, i_remainder_checkersize, mRowCurrent, trace_prep_i : %d, %d, %v, %v\n", i, i%checkersize, mRowCurrent, trace_prep_i) //TODO: why does trace_prep_i change before fill? are ifs executed in parallel? is "=" not one-time?
		if i % checkersize == 0 && i < mColumn-1 { //if i/checkersize remainder is 0, e.g. i=0,3,6,9
			trace_prep_i[i/checkersize] = mRowCurrent //add this mRowCurrent to trace_prep_i before memory is reused, e.g. adding i=0,3,6,9 as trace_prep_i[0,1,2,3][all j]
			fmt.Printf("After fill: i, i_remainder_checkersize, mRowCurrent, trace_prep_i : %d, %d, %v, %v\n", i, i%checkersize, mRowCurrent, trace_prep_i) //TODO
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
		} else if i < mColumn-1 {
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious //reuse mRow variables to save memory, swap slices rather than copy(mRowPrevious,mRowCurrent) for easy update, but only up until the second to last row
		}
	}
	score_highest := mRowCurrent[len(mRowCurrent)-1] //make score_highest variable to be returned
	trace_prep_i[0][0] = 0
	trace_prep_i[0][1] = -430
	trace_prep_i[0][2] = -860
	trace_prep_i[0][3] = -1290
	trace_prep_i[0][4] = -1720 //TODO: remove this chunk after debugging mysterious trace_prep_i fill
	fmt.Printf("score_highest: %d\n", score_highest) //TODO: remove after debugging. PASSED
	fmt.Printf("trace_prep_i: %v\n", trace_prep_i) //TODO: remove after debugging
	fmt.Printf("trace_prep_j: %v\n", trace_prep_j) //TODO: remove after debugging. PASSED

	//start traceback
	var k1, k2, i_inChecker, j_inChecker int //added k1 to cut row into checkersize, and k2 to column into checkersize, internal/within checkerboard i and j, e.g. i_inChecker and j_inChecker can take on values 0,1,2
	route := make([]Cigar, 1)
	for k1, k2 = len(trace_prep_i)-1, len(trace_prep_j)-1; k1 >= 0 && k2 >= 0; { //while loop. TODO: to debug, replace || with &&, so if k1 or k2 reaches 0, no more k loop. Question - Is that the end of traceback?
	//k1, k2 = len(trace_prep_i)-1, len(trace_prep_j)-1 //TODO: try not using for/while loop for k
	//e.g. start from k1=3, while until k1=0, but not at k1-- increment since can skip checkerboards
	//e.g. start from k2=4, go until k2=0, but not at k2-- increment since can skip checkerboards
	//now, fill the prep, aka 4x4 surrounding 3x3 checkerboard
		//mRowPrevious = trace_prep_i[k1] //TODO: uncomment after debugging //initialize the mRowPrevious of the checkerboard. e.g. k1=3, mRowPrevious is now row i=9 (when k1=0, mRowPrevious is row i=0)
		//mRowCurrent[checkersize*k2] = trace_prep_j[k2][checkersize*k1+1] //TODO: uncomment after debugging //e.g. i=k1*checkersize=9, k2=4, mRowCurrent is i=10, mRowCurrent[3*4=12]=trace_prep_j[4][9+1=10] (when k2=0, for each mRowCurrent, will fill j=0 entries)
		//now, fill this checkerboard m and trace from bottom to top, left to right
		if k1==1 && k2==1 { //TODO: remove this chunk after debugging mysterious trace_prep_i and trace_prep_j fill
			mRowPrevious[0] = -1290
			mRowPrevious[1] = -769
			mRowPrevious[2] = -239
			mRowPrevious[3] = 291
			mRowPrevious[4] = -139
			mRowCurrent[3] = -139
		} else if k1==0 && k2==0 { //TODO: remove this chunk after debugging mysterious trace_prep_i and trace_prep_j fill
			mRowPrevious[0] = 0
			mRowPrevious[1] = -430
			mRowPrevious[2] = -860
			mRowPrevious[3] = -1290
			mRowPrevious[4] = -1720
			mRowCurrent[0] = -430
		}
		fmt.Printf("initialize checkerboard: k1, k2, trace_prep_i[k1], mRowPrevious, mRowCurrent: %d, %d, %v, %v, %v\n", k1, k2, trace_prep_i[k1], mRowPrevious, mRowCurrent) //TODO
		for i = checkersize*k1+1; i <= checkersize*(k1+1) && i < mColumn; i++ { //check not only checkersize but also the last chunk <checkersize, break loop if break either condition, e.g. k1=3, start from i=3*3+1=10 since i=3 is prep, need to satisfy i<3*4=12 && i<mColumn=11, so stop at (including) i=10
			i_inChecker = (i-1) % checkersize //e.g. i_inChecker start at 0, and can take on values 0,1,2. Updates with each i
			mRowCurrent[checkersize*k2] = trace_prep_j[k2][checkersize*k1+1+i_inChecker] //this needs to update with i
			for j = checkersize*k2+1; j <= checkersize*(k2+1) && j < len(mRowCurrent); j++ { //e.g. k2=4, start from j=4*3+1=13 since j=12 is prep, need to satisfy j<3*5=15 && j<len(mRowCurrent)=15, so stop at (including) j=14
				j_inChecker = (j-1) % checkersize
				//e.g. i will go from 10-10, j will go from 13-14
				//i==0 and j==0 cases are covered already for mRowCurrent
				if i == 0 && j == 0 {
					trace[i_inChecker][j_inChecker] = 0
				} else if i == 0 {
					trace[i_inChecker][j_inChecker] = 1
				} else if j == 0 {
					trace[i_inChecker][j_inChecker] = 2
				} else {
					mRowCurrent[j], trace[i_inChecker][j_inChecker] = tripleMaxTrace(mRowPrevious[j-1]+scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1]+gapPen, mRowPrevious[j]+gapPen)
					//this should work because mRowPrevious is populated by trace_prep_i, and mRowCurrent[j-1] is populated by trace_prep_j
					//it should be ok even if mRowCurrent isn't completely filled
				}
				fmt.Printf("filled checkerboard: k1, k2, i, j, i_inChecker, j_inChecker, alpha[i-1], beta[j-1], mRowCurrent[j], trace[i_inChecker][j_inChecker]: %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", k1, k2, i, j, i_inChecker, j_inChecker, alpha[i-1], beta[j-1], mRowCurrent[j], trace[i_inChecker][j_inChecker]) //TODO: remove after debugging
				fmt.Printf("mRowPrevious[j-1], scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1], gapPen, mRowPrevious[j]: %d, %d, %d, %d, %d\n", mRowPrevious[j-1], scores[alpha[i-1]][beta[j-1]], mRowCurrent[j-1], gapPen, mRowPrevious[j])
				//TODO: remove fmt after debugging. Magically, now example 1 mCurrent and trace PASSED! Matrix fill success
			}
			if i <= checkersize*(k1+1)-1 && i < mColumn-1 {
				mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
				fmt.Printf("switched\n")
			}
		}

		if (k1 == 0 || k2 == 0) && (i_inChecker == checkersize-1 && j_inChecker == checkersize-1) {
			break
		}

		//TODO: cigar
		//finish traceback and write cigar
		//Consider making this its own function
		var routeIdx int
		for i_inChecker, j_inChecker, routeIdx = checkersize-1, checkersize-1, 0; i_inChecker > 0 && j_inChecker > 0; { //this used to be || but now I think we need &&
			if k1 == len(trace_prep_i)-1 { //if the last bits of the matrix aka biggest k1 and k2, don't access entries that don't exist, start accessing from largest possible i_inChecker and j_inChecker, skipping invalid cases
				i_inChecker = (mColumn-1-1) % checkersize //e.g. largest possible i_inChecker is (11-1-1)%3=0
			}
			if k2 == len(trace_prep_j)-1 {
				j_inChecker = (len(mRowCurrent)-1-1) % checkersize //e.g. largest possible j_inChecker is (15-1-1)%3=1
			}

			if route[routeIdx].RunLength == 0 {
				route[routeIdx].RunLength = 1
				route[routeIdx].Op = trace[i_inChecker][j_inChecker]
			} else if route[routeIdx].Op == trace[i_inChecker][j_inChecker] {
				route[routeIdx].RunLength += 1
			} else {
				route = append(route, Cigar{RunLength: 1, Op: trace[i_inChecker][j_inChecker]})
				routeIdx++
			} //concern: Same Op, Runlength broken up. get rid of route_tmp which was in constGap_lowMem3.go

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
		if i_inChecker == 0 && j_inChecker == 0 { //outside of while loop
			// go to the next checkerboard with lower i, lower j
			k1, k2 = k1-1, k2-1
		} else if i_inChecker == 0 { //but j_inChecker != 0
			// go to the next checkerboard with lower i, same j
			k1 -= 1
		} else if j_inChecker == 0 { //but i_inChecker != 0
			// go to the next checkerboard with lower j, same i
			k2 -= 1
		}
		//Here is where to call the fill matrix function again
	}
	reverseCigar(route)
	//route = append(route, route_tmp...) //use ... to append a []Cigar to another []Cigar

	return score_highest, route
}
