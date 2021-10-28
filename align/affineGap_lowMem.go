package align

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	//"github.com/vertgenlab/gonomics/fasta" //TODO: uncomment later, when use
	"log"
)

// the trace data structure is a 3d slice where the first index is 0,1,2 and represents the match, gap in x (first seq), and gap in y (second seq).
// m used to have the same data structure as trace, but has been simplified into a 2d slice, where the second index for mColumn is removed in order to recycle memory by rows
// for lowMem implementation, have 2 functions initialize Scoring and Trace data structures separately
func initAffineScoring(firstSeqLen int, secondSeqLen int, checkersize_i int, checkersize_j int) ([][]int64, [][]int64, int, [][][]int64, [][][]int64) {
	mRowCurrent := make([][]int64, 3)
	mRowPrevious := make([][]int64, 3)
	var mColumn int = firstSeqLen + 1
	trace_prep_i := make([][][]int64, 3)
	trace_prep_j := make([][][]int64, 3)
	for k := range mRowCurrent { //k ranges through 3 numbers (0,1,2)
		mRowCurrent[k] = make([]int64, secondSeqLen+1)
		mRowPrevious[k] = make([]int64, secondSeqLen+1)
		trace_prep_i[k] = make([][]int64, int(firstSeqLen/checkersize_i)+1) //trace_prep_i saves all rows (i) that are needed to initialize checkerboard tracebacks
		trace_prep_j[k] = make([][]int64, int(secondSeqLen/checkersize_j)+1)  //trace_prep_j saves all columns (j) that are needed to initialize checkerboard tracebacks
		for idx := range trace_prep_i[0] {
			trace_prep_i[k][idx] = make([]int64, secondSeqLen+1)
		}
		for idx := range trace_prep_j[0] {
			trace_prep_j[k][idx] = make([]int64, firstSeqLen+1)
		}
	}
	return mRowCurrent, mRowPrevious, mColumn, trace_prep_i, trace_prep_j
}

func initAffineTrace(firstSeqLen int, secondSeqLen int, checkersize_i int, checkersize_j int) ([][][]ColType) {
	trace_size_i := numbers.Min(firstSeqLen, checkersize_i) //make trace a matrix of size checkersize_i*checkersize_j, unless alpha or beta are shorter, in which case there is no need to allocate a full checkersize of memory
	trace_size_j := numbers.Min(secondSeqLen, checkersize_j)
	trace := make([][][]ColType, 3)
	for k := range trace { //k ranges through 3 numbers (0,1,2)
		trace[k] = make([][]ColType, trace_size_i)
		for i := range trace[0] {
			trace[k][i] = make([]ColType, trace_size_j)
		}
	}
	return trace
}

//like constGap_lowMem, affineGap_lowMem will have 4 steps. Step 1 is different for the variations of affineGap like chunk and multiple, but Step 2-4 can be written together since all variations in affineGap_highMem.go uses the same affineTrace function.
func AffineGap_step1(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64) (int64, int, int, [][][]int64, [][][]int64) {
	//initialize checkersizes here. TODO: make default 10000, and have version where we can change checkersize
	checkersize_i := 3
	checkersize_j := 3

	mRowCurrent, mRowPrevious, mColumn, trace_prep_i, trace_prep_j := initAffineScoring(len(alpha), len(beta), checkersize_i, checkersize_j)

	var i, j int

	for i = 0; i < mColumn; i++ {

		for j = range mRowCurrent[0] {
			if i == 0 && j == 0 {
				mRowCurrent[0][j] = 0
				trace_prep_j[0][j/checkersize_j][i] = mRowCurrent[0][j] //it is implied here that j%checkersize_j==0. It must be saved in trace_prep_j. Although, probably not necessary since initialize all with 0, and affineGap.go doesn't have trace in i==0 && j==0
				mRowCurrent[1][j] = gapOpen
				trace_prep_j[1][j/checkersize_j][i] = mRowCurrent[1][j]
				mRowCurrent[2][j] = gapOpen
				trace_prep_j[2][j/checkersize_j][i] = mRowCurrent[2][j]
			} else if i == 0 {
				mRowCurrent[0][j] = veryNegNum
				mRowCurrent[1][j] = gapExtend + mRowCurrent[1][j-1]
				//trace[1][i][j] = ColI /*new*/
				mRowCurrent[2][j] = veryNegNum
				if j%checkersize_j == 0 {
					trace_prep_j[0][j/checkersize_j][i] = mRowCurrent[0][j]
					trace_prep_j[1][j/checkersize_j][i] = mRowCurrent[1][j]
					trace_prep_j[2][j/checkersize_j][i] = mRowCurrent[2][j]
				}
			} else if j == 0 {
				mRowCurrent[0][j] = veryNegNum
				mRowCurrent[1][j] = veryNegNum
				mRowCurrent[2][j] = gapExtend + mRowPrevious[2][j]
				//trace[2][i][j] = ColD /*new*/
				trace_prep_j[0][j/checkersize_j][i] = mRowCurrent[0][j] //it is implied that j%checkersize_j==0. It must be saved in trace_prep_j
				trace_prep_j[1][j/checkersize_j][i] = mRowCurrent[1][j]
				trace_prep_j[2][j/checkersize_j][i] = mRowCurrent[2][j]
			} else {
				mRowCurrent[0][j], _ = tripleMaxTrace(scores[alpha[i-1]][beta[j-1]]+mRowPrevious[0][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[1][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[2][j-1])
				mRowCurrent[1][j], _ = tripleMaxTrace(gapOpen+gapExtend+mRowCurrent[0][j-1], gapExtend+mRowCurrent[1][j-1], gapOpen+gapExtend+mRowCurrent[2][j-1])
				mRowCurrent[2][j], _ = tripleMaxTrace(gapOpen+gapExtend+mRowPrevious[0][j], gapOpen+gapExtend+mRowPrevious[1][j], gapExtend+mRowPrevious[2][j])
				//mRowCurrent[0][j], trace[0][i][j] = tripleMaxTrace(scores[alpha[i-1]][beta[j-1]]+mRowPrevious[0][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[1][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[2][j-1])
				//mRowCurrent[1][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend+mRowCurrent[0][j-1], gapExtend+mRowCurrent[1][j-1], gapOpen+gapExtend+mRowCurrent[2][j-1])
				//mRowCurrent[2][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend+mRowPrevious[0][j], gapOpen+gapExtend+mRowPrevious[1][j], gapExtend+mRowPrevious[2][j])
				if j%checkersize_j == 0 {
					trace_prep_j[0][j/checkersize_j][i] = mRowCurrent[0][j]
					trace_prep_j[1][j/checkersize_j][i] = mRowCurrent[1][j]
					trace_prep_j[2][j/checkersize_j][i] = mRowCurrent[2][j]
				}
			}
		}

		if i%checkersize_i == 0 && i < mColumn-1 {
			copy(trace_prep_i[0][i/checkersize_i], mRowCurrent[0]) //copy instead of assign values to avoid pointers automatically updating/synchronizing the 2 matrices
			copy(trace_prep_i[1][i/checkersize_i], mRowCurrent[1])
			copy(trace_prep_i[2][i/checkersize_i], mRowCurrent[2])
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
		} else if i < mColumn-1 {
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious //reuse mRow variables to save memory, but only up until the second to last row
		}
	}
	score_highest_i := mColumn - 1
	score_highest_j := len(mRowCurrent[0]) - 1 //TODO: remove this note to self after debugging - the "last J" in affineGap_highMem.go
	score_highest, _ := tripleMaxTrace(mRowCurrent[0][score_highest_j], mRowCurrent[1][score_highest_j], mRowCurrent[2][score_highest_j]) //TODO: remove this note to self after debugging - "maxScore" in affineGap_highMem.go, thd 2nd value tells you direction of "maxScore"
	//TODO: remove these prints
	fmt.Printf("score_highest: %d\n", score_highest)
	fmt.Printf("score_highest_i, score_highest_j: %d, %d\n", score_highest_i, score_highest_j)
	fmt.Printf("trace_prep_i: %v\n", trace_prep_i)
	fmt.Printf("trace_prep_j: %v\n", trace_prep_j)
	return score_highest, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j
}

func AffineGap_step234_testing(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64)([]Cigar) { //input=same as AffineGap_step1 for now, output=route
	//initialize checkersizes here. TODO: make default 10000, and have version where we can change checkersize
	checkersize_i := 3
	checkersize_j := 3

	//Step 1
	//for affineGap, may need a switch here to get the right highest-score output set depending on which variation of affineGap is used
	_, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j := AffineGap_step1(alpha, beta, scores, gapOpen, gapExtend) //the first output score_highest isn't used in this function, so omit

	//Make variables needed for checkerboards
	//k1: the i-index (row-index) of the current checkerboard
	//k2: the j-index (column-index) of the current checkerboard
	//i_inChecker_max: the max i-index (row-index) in the current checkerboard during Step 2 (fillTraceback), either checkersize-1 or a smaller number, e.g. if the alpha/beta sequence to be aligned is not divisible into perfect checkerboards, and the remainder checkerboard is not completely filled
	//j_inChecker_max: the max j-index (column-index) in the current checkerboard during Step 2 (fillTraceback)
	//i_inChecker_min: the min i-index (row-index) in the current checkerboard during Step 3 (writeCigar), either 0 or another number, e.g. if cigar route leaves the checkerboard at a position that is not i_inChecker==0, j_inChecker==0
	//j_inChecker_min: the min j-index (column-index) in the current checkerboard during Step 3 (writeCigar)
	var k1, k2, i_inChecker_max, j_inChecker_max, i_inChecker_min, j_inChecker_min int
	i_inChecker_min = -2                                   //initialize i_inChecker_min != 0, so that the first ever Step 3 (writeCigar) will not interfere with i_inChecker_max
	j_inChecker_min = -2                                   //ditto for j
	trace := initAffineTrace(len(alpha), len(beta), checkersize_i, checkersize_j) //for affineGap, use initAffineTrace function to initialize trace
	route := make([]Cigar, 1) //initialie cigar route and routeIdx
	var routeIdx_current int = 0
	//for affineGap, make a variable to hold k_inChecker_max, outside of for loop
	var k_inChecker_max ColType

	for k1, k2 = int((score_highest_i-1)/checkersize_i), int((score_highest_j-1)/checkersize_j); k1 >= 0 && k2 >= 0; { //use a function of score_highest_i, score_highest_j, and checkersize to initialize the right k1 and k2, go to the correct checkerboard to start traceback

		//Step 2: for a checkerboard, fill traceback, as well as find the max i and j when the checkerboard traceback ends. Since the trace matrix was initialized above, pass it back and forth to recycle the memory instead of creating another trace in Step 2 in every iteration
		//for affineGap, step 2 can be almost the same as constGap, just have 3 dimensions, use gapOpen and gapExtend instead of gapPen, add k_inChecker_max
		trace, i_inChecker_max, j_inChecker_max, k_inChecker_max = fillTraceback_affineGap(alpha, beta, scores, gapOpen, gapExtend, checkersize_i, checkersize_j, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j, k1, k2, i_inChecker_min, j_inChecker_min, trace)

		//Step 3: for a checkerboard, use traceback to write cigar (update route and routeIdx), as well as find the min i and j when the checkerboard cigar ends
		//for affineGap, rewrite step 3 to deal with 3D trace
		route, routeIdx_current, i_inChecker_min, j_inChecker_min = writeCigar_affineGap(trace, i_inChecker_max, j_inChecker_max, route, routeIdx_current, i_inChecker_min, j_inChecker_min, k_inChecker_max)

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
		route, routeIdx_current = lastCigar_affineGap(route, routeIdx_current, 2)
	} else if i_inChecker_min == -1 && j_inChecker_min != -1 { //indicating that Step 3 ended when k1 and k2 reached the smallest combination, and reached i=0, so the last cigar is a "1", e.g. if last cigar entry is the i=0, j=1 square
		route, routeIdx_current = lastCigar_affineGap(route, routeIdx_current, 1)
	} //no more "else" because the only situation left is if Step 3 ended when k1 and k2 reached the smallest combination, and reached both i=0 and j=0, aka the i=0 j=0 square, and there is no sequence there, so no cigar

	//Final processing (reverse route) and return outputs
	reverseCigar(route)
	return route
}

//Step 2
//inputs: the sequences to be aligned alpha and beta, scoring matrix for base matches, penalty for gaps in alignment, checkerboard size for row (i) and column (j), row (i) and column (j) positions of the highest score, trace prep matrices for rows (i) and columns (j)
//				coordinates specifying the current checkerboard in rows (k1) and columns (k2), row (i) and column (j) positions of inChecker_min_Previous which describe where Step 3 (writeCigar) stopped in the previous checkerboard that just went through Step 3 (writeCigar)
//				trace matrix which was initialized in the main ConstGap function and passed back and forth to be recycled
//outputs: trace matrix to be recycled, row (i) and column (j) positions of inChecker_max which describe where Step 2 (fillTraceback) stopped in the current checkerboard that just went through Step 2 (fillTraceback)
//for affineGap, this is 3-dimension version of traceback, return another variable k_inChecker_max_candidate, based on mRowCurrent and j_max
func fillTraceback_affineGap(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64, checkersize_i int, checkersize_j int, score_highest_i int, score_highest_j int, trace_prep_i [][][]int64, trace_prep_j [][][]int64, k1 int, k2 int, i_inChecker_min_Previous int, j_inChecker_min_Previous int, trace [][][]ColType) ([][][]ColType, int, int, ColType) {
	mRowCurrent := make([][]int64, 3)
	mRowPrevious := make([][]int64, 3)
	for k := range mRowCurrent { //k ranges through 3 numbers (0,1,2)
		mRowCurrent[k] = make([]int64, len(beta)+1)
		mRowPrevious[k] = make([]int64, len(beta)+1)
		copy(mRowPrevious[k], trace_prep_i[k][k1])
	}

	//i_inChecker and j_inChecker are internal (within checkerboard) versions of i and j, e.g. for checkersize=3, i_inChecker can take on the values 0, 1, 2
	//i_max and j_max are external (non-inChecker) versions of i_inChecker_max and j_inChecker_max
	var i, j, i_inChecker, j_inChecker, i_inChecker_max, j_inChecker_max, i_max, j_max int
	if i_inChecker_min_Previous >= 0 {
		//use the i_inChecker_min_Previous from the previous k1/k2 loop iteration's Step 3 writeCigar to find the i_max at which to stop filling the trace matrix
		i_max = checkersize_i*k1 + 1 + i_inChecker_min_Previous
	} else {
		//otherwise, since no restriction from previous k1/k2 loop Step 3 writeCigar, find the minimum value between "where the next checkerboard starts" and "the i of the highest score", and that is where to stop filing the trace matrix
		i_max = numbers.Min(checkersize_i*(k1+1), score_highest_i)
	}
	if j_inChecker_min_Previous >= 0 {
		j_max = checkersize_j*k2 + 1 + j_inChecker_min_Previous
	} else {
		j_max = numbers.Min(checkersize_j*(k2+1), score_highest_j)
	}
	i_inChecker_max = (i_max - 1) % checkersize_i
	j_inChecker_max = (j_max - 1) % checkersize_j

	for i = checkersize_i*k1 + 1; i <= i_max; i++ {
		i_inChecker = (i - 1) % checkersize_i
		mRowCurrent[0][checkersize_j*k2] = trace_prep_j[0][k2][checkersize_j*k1+1+i_inChecker]
		mRowCurrent[1][checkersize_j*k2] = trace_prep_j[1][k2][checkersize_j*k1+1+i_inChecker]
		mRowCurrent[2][checkersize_j*k2] = trace_prep_j[2][k2][checkersize_j*k1+1+i_inChecker]

		for j = checkersize_j*k2 + 1; j <= j_max; j++ {
			j_inChecker = (j - 1) % checkersize_j
			mRowCurrent[0][j], trace[0][i_inChecker][j_inChecker] = tripleMaxTrace(scores[alpha[i-1]][beta[j-1]]+mRowPrevious[0][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[1][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[2][j-1])
			mRowCurrent[1][j], trace[1][i_inChecker][j_inChecker] = tripleMaxTrace(gapOpen+gapExtend+mRowCurrent[0][j-1], gapExtend+mRowCurrent[1][j-1], gapOpen+gapExtend+mRowCurrent[2][j-1])
			mRowCurrent[2][j], trace[2][i_inChecker][j_inChecker] = tripleMaxTrace(gapOpen+gapExtend+mRowPrevious[0][j], gapOpen+gapExtend+mRowPrevious[1][j], gapExtend+mRowPrevious[2][j])
			//it is ok even if mRowCurrent isn't completely filled, aka only part of mRowCurrent is needed for the checkerboard traceback and writing cigar
		}

		if i <= checkersize_i*(k1+1)-1 && i <= score_highest_i-1 {
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
		}
	}

	//for affineGap, return another variable k_inChecker_max, based on mRowCurrent and j_max
	_, k_inChecker_max := tripleMaxTrace(mRowCurrent[0][j_max], mRowCurrent[1][j_max], mRowCurrent[2][j_max])

	return trace, i_inChecker_max, j_inChecker_max, k_inChecker_max
}

//Step 3
//inputs: trace matrix, row (i) and column (j) positions of inChecker_max_Previous which describe where Step 2 (fillTraceback) stopped in the previous checkerboard that just went through Step 2 (fillTraceback)
//				the growing route describing the collection of cigars of the entire alignment, the index of the cigar that is currently being built, row (i) and column (j) positions of inChecker_min_Previous which describe where Step 3 (writeCigar) stopped in the previous checkerboard that just went through Step 3 (writeCigar)
//outputs: the updated grown route describine the collection of cigars of the entire alignment, the updated index of the cigar that is currently being built, row (i) and column (j) positions of inChecker_min which describe where Step 3 (writeCigar) stopped in the current checkerboard that just went through Step 3 (writeCigar)
//for affineGap, made changes for 3D trace
func writeCigar_affineGap(trace [][][]ColType, i_inChecker_max_Previous int, j_inChecker_max_Previous int, route []Cigar, routeIdx_current int, i_inChecker_min_Previous int, j_inChecker_min_Previous int, k_inChecker_max ColType) ([]Cigar, int, int, int) {

	var i_inChecker, j_inChecker, routeIdx, i_inChecker_min, j_inChecker_min, i_inChecker_max, j_inChecker_max int
	route_updated := route

	//use the i_inChecker_min_Previous from the previous k1/k2 loop iteration's Step 3 writeCigar to find the i_inChecker_max at which to start writing cigar
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
			route_updated[routeIdx].Op = k_inChecker_max //trace[i_inChecker][j_inChecker] in constGap.go
		} else if route_updated[routeIdx].Op == k_inChecker_max { //trace[i_inChecker][j_inChecker] in constGap.go
			route_updated[routeIdx].RunLength += 1
		} else {
			route_updated = append(route_updated, Cigar{RunLength: 1, Op: k_inChecker_max}) //trace[i_inChecker][j_inChecker] in constGap.go
			routeIdx++
		}

		//update i_inChecker and j_inChecker
		switch k_inChecker_max { //trace[i_inChecker][j_inChecker] in constGap.go
		case 0:
			k_inChecker_max = trace[k_inChecker_max][i_inChecker][j_inChecker] //k = trace[k][i][j] in affineGap.go. This is a step to update k_inChecker_max unique to affineGap.go
			i_inChecker, j_inChecker = i_inChecker-1, j_inChecker-1
		case 1:
			k_inChecker_max = trace[k_inChecker_max][i_inChecker][j_inChecker]
			j_inChecker -= 1
		case 2:
			k_inChecker_max = trace[k_inChecker_max][i_inChecker][j_inChecker]
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
//inputs: the growing route describing the collection of cigars of the entire alignment except the last cigar, the index of the cigar that is currently being built, the ColType (M=0,I=1,D=2) of the last cigar entry
//outputs: the updated final route describing the collection of cigars of the entire alignment including the last cigar, the updated and index of the cigar that is currently being built
func lastCigar_affineGap(route []Cigar, routeIdx_current int, Op_end ColType) ([]Cigar, int) {

	if route[routeIdx_current].Op == Op_end {
		route[routeIdx_current].RunLength += 1
	} else {
		route = append(route, Cigar{RunLength: 1, Op: Op_end})
		routeIdx_current++
	}

	return route, routeIdx_current
}
