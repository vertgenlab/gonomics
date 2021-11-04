package align

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

//reminder that unlike constGap, affineGap is an alignment algorithm that favors having a long gap compared to multiple gaps. In the code, affineGap data structures have an additional "k" dimension, which can take on values 0,1,2 and represent the "alternative universes" of scores where the correct alignment produces Match(0), Insertion(1) and Deletion(2)

//"Step 0"
//make data structures for scoring
// the trace data structure is a 3d slice where the first index is 0,1,2 and represents the match, gap in x (first seq), and gap in y (second seq).
// m used to have the same data structure as trace, but has been simplified into a 2d slice, where the second index for mColumn is removed in order to recycle memory by rows
// for lowMem4 checkerboard implementation, have 2 functions initialize Scoring and Trace data structures separately
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

//"Step 0"
//make data structures for tracing
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

//This version of AffineGap has a fixed checkersize of 10000*10000
func AffineGap(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64)(int64, []Cigar) {
	var checkersize_i, checkersize_j int
	checkersize_i = 10000
	checkersize_j = 10000
	
	//Step 1: find highest score, as well as get the position (i and j) of the highest score, and materials needed to fill traceback and write cigar in checkerboards
	score_highest, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j := highestScore_affineGap(alpha, beta, scores, gapOpen, gapExtend, checkersize_i, checkersize_j)

	//Make variables needed for checkerboards
	//k1: the i-index (row-index) of the current checkerboard
	//k2: the j-index (column-index) of the current checkerboard
	//i_inChecker_max: the max i-index (row-index) in the current checkerboard during Step 2 (fillTraceback), either checkersize-1 or a smaller number, e.g. if the alpha/beta sequence to be aligned is not divisible into perfect checkerboards, and the remainder checkerboard is not completely filled
	//j_inChecker_max: the max j-index (column-index) in the current checkerboard during Step 2 (fillTraceback)
	//k_inChecker_max: the k-index (out of 0,1,2) of the max i and j in the current checkerboard during Step 2 (fillTraceback)
	//i_inChecker_min: the min i-index (row-index) in the current checkerboard during Step 3 (writeCigar), either 0 or another number, e.g. if cigar route leaves the checkerboard at a position that is not i_inChecker==0, j_inChecker==0
	//j_inChecker_min: the min j-index (column-index) in the current checkerboard during Step 3 (writeCigar)
	//k_inChecker_min: the k-index (out of 0,1,2) of the min i and j in the current checkerboard during Step 3 (writeCigar)
	//note that k1 amd k2 (coordinates for checkerboads) are not related to k_inChecker_max and k_inChecker_min (out of 0,1,2). Also note that k_inCkecer_max and k_inChecker_min do not determine traceback path (e.g. determine the next k1 and k2 to go to), but influence cigar route
	var k1, k2, i_inChecker_max, j_inChecker_max, i_inChecker_min, j_inChecker_min int
	i_inChecker_min = -2 //initialize i_inChecker_min != 0, so that the first ever Step 3 (writeCigar) will not interfere with i_inChecker_max
	j_inChecker_min = -2 //ditto for j
	trace := initAffineTrace(len(alpha), len(beta), checkersize_i, checkersize_j) //for affineGap, use initAffineTrace function to initialize trace
	route := make([]Cigar, 1) //initialie cigar route and routeIdx
	var routeIdx_current int = 0
	var k_inChecker_max, k_inChecker_min ColType //for affineGap, make a variable to hold k_inChecker_max, outside of for loop

	for k1, k2 = int((score_highest_i-1)/checkersize_i), int((score_highest_j-1)/checkersize_j); k1 >= 0 && k2 >= 0; { //use a function of score_highest_i, score_highest_j, and checkersize to initialize the right k1 and k2, go to the correct checkerboard to start traceback

		//Step 2: for a checkerboard, fill traceback, as well as find the max i and j when the checkerboard traceback ends. Since the trace matrix was initialized above, pass it back and forth to recycle the memory instead of creating another trace in Step 2 in every iteration
		//for affineGap, step 2 is almost the same as constGap, just have 3 dimensions of k (0,1,2), use gapOpen and gapExtend instead of gapPen, add k_inChecker_max
		trace, i_inChecker_max, j_inChecker_max, k_inChecker_max = fillTraceback_affineGap(alpha, beta, scores, gapOpen, gapExtend, checkersize_i, checkersize_j, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j, k1, k2, i_inChecker_min, j_inChecker_min, trace)

		//Step 3: for a checkerboard, use traceback to write cigar (update route and routeIdx), as well as find the min i and j when the checkerboard cigar ends
		//for affineGap, step 3 is almost the same as constGap, just have 3 dimensions of k (0,1,2), add k_inChecker_min
		route, routeIdx_current, i_inChecker_min, j_inChecker_min, k_inChecker_min = writeCigar_affineGap(trace, i_inChecker_max, j_inChecker_max, k_inChecker_max, route, routeIdx_current, i_inChecker_min, j_inChecker_min, k_inChecker_min)

		//Use Step 3's i_inChecker_min and j_inChecker_min to find the next checkerboard to go to and where to start in that checkerboard (update k1, k2)
		if i_inChecker_min < 0 && j_inChecker_min < 0 {
			k1, k2 = k1-1, k2-1 //go to the next checkerboard with lower i, lower j
		} else if i_inChecker_min < 0 { //but j_inChecker != 0
			k1 -= 1 //go to the next checkerboard with lower i, same j
		} else if j_inChecker_min < 0 { //but i_inChecker != 0
			k2 -= 1 //go to the next checkerboard with lower j, same i
		}

	}

	//Step 4: write the last cigar entry. Since this step doesn't need special features for affineGap, reuse constGap function
	//This step is necessary because the row i=0 and the column j=0 are always stored in trace_prep, and never filled in in any checkerboard, but the cigar ends by the route going into 1 box/entry in either the row i=0 or the column j=0
	if i_inChecker_min != -1 && j_inChecker_min == -1 { //indicating that Step 3 ended when k1 and k2 reached the smallest combination, and reached j=0, so the last cigar is a "2", e.g. if last cigar entry is the i=1, j=0 square
		route, routeIdx_current = lastCigar(route, routeIdx_current, 2)
	} else if i_inChecker_min == -1 && j_inChecker_min != -1 { //indicating that Step 3 ended when k1 and k2 reached the smallest combination, and reached i=0, so the last cigar is a "1", e.g. if last cigar entry is the i=0, j=1 square
		route, routeIdx_current = lastCigar(route, routeIdx_current, 1)
	} //no more "else" because the only situation left is if Step 3 ended when k1 and k2 reached the smallest combination, and reached both i=0 and j=0, aka the i=0 j=0 square, and there is no sequence there, so no cigar

	//Final processing (reverse route) and return outputs
	reverseCigar(route)
	return score_highest, route
}

//This version of AffineGap needs additional inputs and allows customization of checkersize_i and checkersize_j
func AffineGap_customizeCheckersize(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64, checkersize_i int, checkersize_j int)(int64, []Cigar) { //input=same as AffineGap_step1 for now, output=route
	//Step 1: find highest score, as well as get the position (i and j) of the highest score, and materials needed to fill traceback and write cigar in checkerboards
	score_highest, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j := highestScore_affineGap(alpha, beta, scores, gapOpen, gapExtend, checkersize_i, checkersize_j)

	//Make variables needed for checkerboards
	//k1: the i-index (row-index) of the current checkerboard
	//k2: the j-index (column-index) of the current checkerboard
	//i_inChecker_max: the max i-index (row-index) in the current checkerboard during Step 2 (fillTraceback), either checkersize-1 or a smaller number, e.g. if the alpha/beta sequence to be aligned is not divisible into perfect checkerboards, and the remainder checkerboard is not completely filled
	//j_inChecker_max: the max j-index (column-index) in the current checkerboard during Step 2 (fillTraceback)
	//k_inChecker_max: the k-index (out of 0,1,2) of the max i and j in the current checkerboard during Step 2 (fillTraceback)
	//i_inChecker_min: the min i-index (row-index) in the current checkerboard during Step 3 (writeCigar), either 0 or another number, e.g. if cigar route leaves the checkerboard at a position that is not i_inChecker==0, j_inChecker==0
	//j_inChecker_min: the min j-index (column-index) in the current checkerboard during Step 3 (writeCigar)
	//k_inChecker_min: the k-index (out of 0,1,2) of the min i and j in the current checkerboard during Step 3 (writeCigar)
	//note that k1 amd k2 (coordinates for checkerboads) are not related to k_inChecker_max and k_inChecker_min (out of 0,1,2). Also note that k_inCkecer_max and k_inChecker_min do not determine traceback path (e.g. determine the next k1 and k2 to go to), but influence cigar route
	var k1, k2, i_inChecker_max, j_inChecker_max, i_inChecker_min, j_inChecker_min int
	i_inChecker_min = -2 //initialize i_inChecker_min != 0, so that the first ever Step 3 (writeCigar) will not interfere with i_inChecker_max
	j_inChecker_min = -2 //ditto for j
	trace := initAffineTrace(len(alpha), len(beta), checkersize_i, checkersize_j) //for affineGap, use initAffineTrace function to initialize trace
	route := make([]Cigar, 1) //initialie cigar route and routeIdx
	var routeIdx_current int = 0
	var k_inChecker_max, k_inChecker_min ColType //for affineGap, make a variable to hold k_inChecker_max, outside of for loop

	for k1, k2 = int((score_highest_i-1)/checkersize_i), int((score_highest_j-1)/checkersize_j); k1 >= 0 && k2 >= 0; { //use a function of score_highest_i, score_highest_j, and checkersize to initialize the right k1 and k2, go to the correct checkerboard to start traceback

		//Step 2: for a checkerboard, fill traceback, as well as find the max i and j when the checkerboard traceback ends. Since the trace matrix was initialized above, pass it back and forth to recycle the memory instead of creating another trace in Step 2 in every iteration
		//for affineGap, step 2 is almost the same as constGap, just have 3 dimensions of k (0,1,2), use gapOpen and gapExtend instead of gapPen, add k_inChecker_max
		trace, i_inChecker_max, j_inChecker_max, k_inChecker_max = fillTraceback_affineGap(alpha, beta, scores, gapOpen, gapExtend, checkersize_i, checkersize_j, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j, k1, k2, i_inChecker_min, j_inChecker_min, trace)

		//Step 3: for a checkerboard, use traceback to write cigar (update route and routeIdx), as well as find the min i and j when the checkerboard cigar ends
		//for affineGap, step 3 is almost the same as constGap, just have 3 dimensions of k (0,1,2), add k_inChecker_min
		route, routeIdx_current, i_inChecker_min, j_inChecker_min, k_inChecker_min = writeCigar_affineGap(trace, i_inChecker_max, j_inChecker_max, k_inChecker_max, route, routeIdx_current, i_inChecker_min, j_inChecker_min, k_inChecker_min)

		//Use Step 3's i_inChecker_min and j_inChecker_min to find the next checkerboard to go to and where to start in that checkerboard (update k1, k2)
		if i_inChecker_min < 0 && j_inChecker_min < 0 {
			k1, k2 = k1-1, k2-1 //go to the next checkerboard with lower i, lower j
		} else if i_inChecker_min < 0 { //but j_inChecker != 0
			k1 -= 1 //go to the next checkerboard with lower i, same j
		} else if j_inChecker_min < 0 { //but i_inChecker != 0
			k2 -= 1 //go to the next checkerboard with lower j, same i
		}

	}

	//Step 4: write the last cigar entry. Since this step doesn't need special features for affineGap, reuse constGap function
	//This step is necessary because the row i=0 and the column j=0 are always stored in trace_prep, and never filled in in any checkerboard, but the cigar ends by the route going into 1 box/entry in either the row i=0 or the column j=0
	if i_inChecker_min != -1 && j_inChecker_min == -1 { //indicating that Step 3 ended when k1 and k2 reached the smallest combination, and reached j=0, so the last cigar is a "2", e.g. if last cigar entry is the i=1, j=0 square
		route, routeIdx_current = lastCigar(route, routeIdx_current, 2)
	} else if i_inChecker_min == -1 && j_inChecker_min != -1 { //indicating that Step 3 ended when k1 and k2 reached the smallest combination, and reached i=0, so the last cigar is a "1", e.g. if last cigar entry is the i=0, j=1 square
		route, routeIdx_current = lastCigar(route, routeIdx_current, 1)
	} //no more "else" because the only situation left is if Step 3 ended when k1 and k2 reached the smallest combination, and reached both i=0 and j=0, aka the i=0 j=0 square, and there is no sequence there, so no cigar

	//Final processing (reverse route) and return outputs
	reverseCigar(route)
	return score_highest, route
}

//Step 1
//inputs: the sequences to be aligned alpha and beta, scoring matrix for base matches, penalty for opening and extending gaps in alignment, checkerboard size for row (i) and column (j)
//outputs: highest score from the alignment, row (i) and column (j) positions of the highest score, trace_prep matrices for rows (i) and columns (j)
func highestScore_affineGap(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64, checkersize_i int, checkersize_j int) (int64, int, int, [][][]int64, [][][]int64) {
	//initialize data structures for getting highest score
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
				trace_prep_j[0][j/checkersize_j][i] = mRowCurrent[0][j] //it is implied that j%checkersize_j==0. It must be saved in trace_prep_j
				trace_prep_j[1][j/checkersize_j][i] = mRowCurrent[1][j]
				trace_prep_j[2][j/checkersize_j][i] = mRowCurrent[2][j]
			} else {
				mRowCurrent[0][j], _ = tripleMaxTrace(scores[alpha[i-1]][beta[j-1]]+mRowPrevious[0][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[1][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[2][j-1])
				mRowCurrent[1][j], _ = tripleMaxTrace(gapOpen+gapExtend+mRowCurrent[0][j-1], gapExtend+mRowCurrent[1][j-1], gapOpen+gapExtend+mRowCurrent[2][j-1])
				mRowCurrent[2][j], _ = tripleMaxTrace(gapOpen+gapExtend+mRowPrevious[0][j], gapOpen+gapExtend+mRowPrevious[1][j], gapExtend+mRowPrevious[2][j])
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
	score_highest_j := len(mRowCurrent[0]) - 1 //score_highest_j was denoted with "last J" in affineGap_highMem.go
	score_highest, _ := tripleMaxTrace(mRowCurrent[0][score_highest_j], mRowCurrent[1][score_highest_j], mRowCurrent[2][score_highest_j]) //denoted by "maxScore" in affineGap_highMem.go. The 2nd output variable is its direction (0,1,2)
	return score_highest, score_highest_i, score_highest_j, trace_prep_i, trace_prep_j
}

//Step 2
//inputs: the sequences to be aligned alpha and beta, scoring matrix for base matches, penalty for opening and extending gaps in alignment, checkerboard size for row (i) and column (j), row (i) and column (j) positions of the highest score, trace prep matrices for rows (i) and columns (j)
//				coordinates specifying the current checkerboard in rows (k1) and columns (k2), row (i) and column (j) positions of inChecker_min_Previous which describe where Step 3 (writeCigar) stopped in the previous checkerboard that just went through Step 3 (writeCigar)
//				trace matrix which was initialized in the main ConstGap function and passed back and forth to be recycled
//outputs: trace matrix to be recycled, row (i), column (j) and dimension (k) positions of inChecker_max which describe where Step 2 (fillTraceback) stopped in the current checkerboard that just went through Step 2 (fillTraceback)
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
//inputs: trace matrix, row (i), column (j) and dimension (k) positions of inChecker_max_Previous which describe where Step 2 (fillTraceback) stopped in the previous checkerboard that just went through Step 2 (fillTraceback)
//				the growing route describing the collection of cigars of the entire alignment, the index of the cigar that is currently being built, row (i), column (j) and dimension (k) positions of inChecker_min_Previous which describe where Step 3 (writeCigar) stopped in the previous checkerboard that just went through Step 3 (writeCigar)
//outputs: the updated grown route describine the collection of cigars of the entire alignment, the updated index of the cigar that is currently being built, row (i), column (j) and dimension (k) positions of inChecker_min which describe where Step 3 (writeCigar) stopped in the current checkerboard that just went through Step 3 (writeCigar)
func writeCigar_affineGap(trace [][][]ColType, i_inChecker_max_Previous int, j_inChecker_max_Previous int, k_inChecker_max_Previous ColType, route []Cigar, routeIdx_current int, i_inChecker_min_Previous int, j_inChecker_min_Previous int, k_inChecker_min_Previous ColType) ([]Cigar, int, int, int, ColType) {

	var i_inChecker, j_inChecker, routeIdx, i_inChecker_min, j_inChecker_min, i_inChecker_max, j_inChecker_max int
	var k_inChecker, k_inChecker_min ColType

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
	if i_inChecker_min_Previous >= 0 && j_inChecker_max_Previous >= 0 {
		k_inChecker = k_inChecker_min_Previous
	} else {
		k_inChecker = k_inChecker_max_Previous
	}

	for i_inChecker, j_inChecker, routeIdx = i_inChecker_max, j_inChecker_max, routeIdx_current; i_inChecker >= 0 && j_inChecker >= 0; {

		//write cigar segment
		if route_updated[routeIdx].RunLength == 0 {
			route_updated[routeIdx].RunLength = 1
			route_updated[routeIdx].Op = k_inChecker
		} else if route_updated[routeIdx].Op == k_inChecker {
			route_updated[routeIdx].RunLength += 1
		} else {
			route_updated = append(route_updated, Cigar{RunLength: 1, Op: k_inChecker})
			routeIdx++
		}

		//update i_inChecker and j_inChecker
		switch k_inChecker {
		case 0:
			k_inChecker = trace[k_inChecker][i_inChecker][j_inChecker]
			i_inChecker, j_inChecker = i_inChecker-1, j_inChecker-1
		case 1:
			k_inChecker = trace[k_inChecker][i_inChecker][j_inChecker]
			j_inChecker -= 1
		case 2:
			k_inChecker = trace[k_inChecker][i_inChecker][j_inChecker]
			i_inChecker -= 1
		default:
			log.Fatalf("Error: unexpected traceback")
		}

		i_inChecker_min = i_inChecker
		j_inChecker_min = j_inChecker
		k_inChecker_min = k_inChecker

	}

	return route_updated, routeIdx, i_inChecker_min, j_inChecker_min, k_inChecker_min //return the routeIdx that was updated in the loop, not the function input "routeIdx_current"
}
