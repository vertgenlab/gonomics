package align

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	//"github.com/vertgenlab/gonomics/fasta" //TODO: uncomment later, when use
	//"log" //TODO: uncomment later, when use
)

// the trace data structure is a 3d slice where the first index is 0,1,2 and represents the match, gap in x (first seq), and gap in y (second seq).
// m used to have the same data structure as trace, but has been simplified into a 2d slice, where the second index for mColumn is removed in order to recycle memory by rows
// for lowMem implementation, have 2 functions initialize Scoring and Trace data structures separately
func initAffineScoringAndTrace_testing(firstSeqLen int, secondSeqLen int, checkersize_i int, checkersize_j int) ([][]int64, [][]int64, int, [][][]int64, [][][]int64) {
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

func initAffineTrace_testing(firstSeqLen int, secondSeqLen int, checkersize_i int, checkersize_j int) ([][][]ColType) {
		trace := make([][][]ColType, 3)
		for k := range trace { //k ranges through 3 numbers (0,1,2)
			trace[k] = make([][]ColType, firstSeqLen+1)
			for i := range trace[0] {
				trace[k][i] = make([]ColType, secondSeqLen+1)
			}
		}
		return trace
}

//like constGap_lowMem, affineGap_lowMem will have 4 steps. Step 1 is different for the variations of affineGap like chunk and multiple, but Step 2-4 can be written together since all variations in affineGap_highMem.go uses the same affineTrace function.
func AffineGap_step1_testing(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64) (int64) { //TODO: replace output with actual output data types (int64, int, int, ColType, [][][]int64, [][][]int64)
	//initialize checkersizes here. TODO: make default 10000, and have version where we can change checkersize
	checkersize_i := 3
	checkersize_j := 3

	mRowCurrent, mRowPrevious, mColumn, trace_prep_i, trace_prep_j := initAffineScoringAndTrace_testing(len(alpha), len(beta), checkersize_i, checkersize_j)

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
	score_highest, score_highest_k := tripleMaxTrace(mRowCurrent[0][score_highest_j], mRowCurrent[1][score_highest_j], mRowCurrent[2][score_highest_j]) //TODO: remove this note to self after debugging - "maxScore" in affineGap_highMem.go, tell you direction of "maxScore" in score_highest_k

	//TODO: make these prints become output
	fmt.Printf("score_highest: %d\n", score_highest)
	fmt.Printf("score_highest_i, score_highest_j, score_highest_k: %d, %d, %d\n", score_highest_i, score_highest_j, score_highest_k)
	fmt.Printf("trace_prep_i: %v\n", trace_prep_i)
	fmt.Printf("trace_prep_j: %v\n", trace_prep_j)
	return score_highest
}
