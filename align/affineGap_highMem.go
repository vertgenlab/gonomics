package align

import (
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

// the trace data structure is a 3d slice where the first index is 0,1,2 and represents the match, gap in x (first seq), and gap in y (second seq).
// m used to have the same data structure as trace, but has been simplified into a 2d slice, where the second index for mColumn is removed in order to recycle memory by rows.
func initAffineScoringAndTrace(firstSeqLen int, secondSeqLen int) ([][]int64, [][]int64, int, [][][]ColType) {
	mRowCurrent := make([][]int64, 3)
	mRowPrevious := make([][]int64, 3)
	var mColumn int = firstSeqLen + 1
	trace := make([][][]ColType, 3)
	for k := range trace { //k ranges through 3 numbers (0,1,2)
		trace[k] = make([][]ColType, firstSeqLen+1)
		mRowCurrent[k] = make([]int64, secondSeqLen+1)
		mRowPrevious[k] = make([]int64, secondSeqLen+1)
		for i := range trace[0] {
			trace[k][i] = make([]ColType, secondSeqLen+1)
		}
	}
	return mRowCurrent, mRowPrevious, mColumn, trace
}

func initAffineScoringAndTraceRecycle(firstSeqLen int, secondSeqLen int, mRowCurrent, mRowPrevious [][]int64, trace [][][]ColType) ([][]int64, [][]int64, int, [][][]ColType) {
	mRowCurrent = setSliceSize(mRowCurrent, 3)
	mRowPrevious = setSliceSize(mRowPrevious, 3)
	var mColumn int = firstSeqLen + 1
	trace = setSliceSize(trace, 3)
	for k := range trace { //k ranges through 3 numbers (0,1,2)
		trace[k] = setSliceSize(trace[k], firstSeqLen+1)
		mRowCurrent[k] = setSliceSize(mRowCurrent[k], secondSeqLen+1)
		mRowPrevious[k] = setSliceSize(mRowPrevious[k], secondSeqLen+1)
		for i := range trace[0] {
			trace[k][i] = setSliceSize(trace[k][i], secondSeqLen+1)
		}
	}
	return mRowCurrent, mRowPrevious, mColumn, trace
}

func setSliceSize[t any](s []t, size int) []t {
	if len(s) == size {
		return s
	}

	if cap(s) >= size {
		return s[:size]
	}

	return make([]t, size)
}

func affineTrace(mRowCurrent [][]int64, mColumn int, trace [][][]ColType) (int64, []Cigar) {
	route := make([]Cigar, 1)
	lastI := mColumn - 1             //the last I
	lastJ := len(mRowCurrent[0]) - 1 //the last J
	maxScore, k := tripleMaxTrace(mRowCurrent[0][lastJ], mRowCurrent[1][lastJ], mRowCurrent[2][lastJ])
	for i, j, routeIdx := lastI, lastJ, 0; i > 0 || j > 0; {
		if route[routeIdx].RunLength == 0 {
			route[routeIdx].RunLength = 1
			route[routeIdx].Op = k
		} else if route[routeIdx].Op == k {
			route[routeIdx].RunLength += 1
		} else {
			route = append(route, Cigar{RunLength: 1, Op: k})
			routeIdx++
		}
		switch k {
		case ColM:
			k = trace[k][i][j]
			i--
			j--
		case ColI:
			k = trace[k][i][j]
			j--
		case ColD:
			k = trace[k][i][j]
			i--
		default:
			log.Fatalf("Error: unexpected traceback")
		}
	}
	reverseCigar(route)
	return maxScore, route
}

func expandCigarRunLength(route []Cigar, chunkSize int64) {
	for i := range route {
		route[i].RunLength *= chunkSize
	}
}

func AffineGap_highMem(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64) (int64, []Cigar) {
	return affineGap_highMem(alpha, beta, scores, gapOpen, gapExtend, false)
}

// AffineGapLocal functions identically to AffineGap_highMem, but it does not penalize for gaps placed at the beginning or end of the alignment.
// This property enables AffineGap to be used for local alignment, such as aligning a 150bp sequencing read to a 1kb reference sequence.
func AffineGapLocal(target []dna.Base, query []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64) (int64, []Cigar) {
	return affineGap_highMem(target, query, scores, gapOpen, gapExtend, true)
}

type TargetQueryPair struct {
	Target []dna.Base
	Query  []dna.Base
	Score  int64
	Cigar  []Cigar
}

func GoAffineGapLocalEngine(scores [][]int64, gapOpen int64, gapExtend int64) (inputs chan<- TargetQueryPair, outputs <-chan TargetQueryPair) {
	i := make(chan TargetQueryPair, 1000)
	o := make(chan TargetQueryPair, 1000)
	go goAffineGap_highMem(i, o, scores, gapOpen, gapExtend, true)
	return i, o
}

func goAffineGap_highMem(inputs <-chan TargetQueryPair, outputs chan<- TargetQueryPair, scores [][]int64, gapOpen int64, gapExtend int64, freeEndGaps bool) {
	var mRowCurrent, mRowPrevious [][]int64
	var mColumn int
	var trace [][][]ColType
	var alpha, beta []dna.Base
	for in := range inputs {
		alpha = in.Target
		beta = in.Query
		mRowCurrent, mRowPrevious, mColumn, trace = initAffineScoringAndTraceRecycle(len(alpha), len(beta), mRowCurrent, mRowPrevious, trace)
		for i := 0; i < mColumn; i++ {
			for j := range mRowCurrent[0] {
				if i == 0 && j == 0 {
					mRowCurrent[0][j] = 0
					mRowCurrent[1][j] = gapOpen
					if freeEndGaps {
						mRowCurrent[2][j] = 0
					} else {
						mRowCurrent[2][j] = gapOpen
					}
				} else if i == 0 {
					mRowCurrent[0][j] = veryNegNum
					mRowCurrent[1][j] = gapExtend + mRowCurrent[1][j-1]
					trace[1][i][j] = ColI
					mRowCurrent[2][j] = veryNegNum
				} else if j == 0 {
					mRowCurrent[0][j] = veryNegNum
					mRowCurrent[1][j] = veryNegNum
					if freeEndGaps {
						mRowCurrent[2][j] = 0 + mRowPrevious[2][j]
					} else {
						mRowCurrent[2][j] = gapExtend + mRowPrevious[2][j]
					}
					trace[2][i][j] = ColD
				} else if freeEndGaps && j == len(mRowCurrent[0])-1 {
					mRowCurrent[0][j], trace[0][i][j] = tripleMaxTrace(scores[alpha[i-1]][beta[j-1]]+mRowPrevious[0][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[1][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[2][j-1])
					mRowCurrent[1][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend+mRowCurrent[0][j-1], gapExtend+mRowCurrent[1][j-1], gapOpen+gapExtend+mRowCurrent[2][j-1])
					mRowCurrent[2][j], trace[2][i][j] = tripleMaxTrace(0+0+mRowPrevious[0][j], 0+0+mRowPrevious[1][j], 0+mRowPrevious[2][j])
				} else {
					mRowCurrent[0][j], trace[0][i][j] = tripleMaxTrace(scores[alpha[i-1]][beta[j-1]]+mRowPrevious[0][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[1][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[2][j-1])
					mRowCurrent[1][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend+mRowCurrent[0][j-1], gapExtend+mRowCurrent[1][j-1], gapOpen+gapExtend+mRowCurrent[2][j-1])
					mRowCurrent[2][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend+mRowPrevious[0][j], gapOpen+gapExtend+mRowPrevious[1][j], gapExtend+mRowPrevious[2][j])
				}
			}
			if i < mColumn-1 {
				mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
			}
		}
		maxScore, route := affineTrace(mRowCurrent, mColumn, trace)
		in.Score, in.Cigar = maxScore, route
		outputs <- in
	}
	close(outputs)
}

func affineGap_highMem(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64, freeEndGaps bool) (int64, []Cigar) {
	mRowCurrent, mRowPrevious, mColumn, trace := initAffineScoringAndTrace(len(alpha), len(beta))
	for i := 0; i < mColumn; i++ {
		for j := range mRowCurrent[0] {
			if i == 0 && j == 0 {
				mRowCurrent[0][j] = 0
				mRowCurrent[1][j] = gapOpen
				if freeEndGaps {
					mRowCurrent[2][j] = 0
				} else {
					mRowCurrent[2][j] = gapOpen
				}
			} else if i == 0 {
				mRowCurrent[0][j] = veryNegNum
				mRowCurrent[1][j] = gapExtend + mRowCurrent[1][j-1]
				trace[1][i][j] = ColI
				mRowCurrent[2][j] = veryNegNum
			} else if j == 0 {
				mRowCurrent[0][j] = veryNegNum
				mRowCurrent[1][j] = veryNegNum
				if freeEndGaps {
					mRowCurrent[2][j] = 0 + mRowPrevious[2][j]
				} else {
					mRowCurrent[2][j] = gapExtend + mRowPrevious[2][j]
				}
				trace[2][i][j] = ColD
			} else if freeEndGaps && j == len(mRowCurrent[0])-1 {
				mRowCurrent[0][j], trace[0][i][j] = tripleMaxTrace(scores[alpha[i-1]][beta[j-1]]+mRowPrevious[0][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[1][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[2][j-1])
				mRowCurrent[1][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend+mRowCurrent[0][j-1], gapExtend+mRowCurrent[1][j-1], gapOpen+gapExtend+mRowCurrent[2][j-1])
				mRowCurrent[2][j], trace[2][i][j] = tripleMaxTrace(0+0+mRowPrevious[0][j], 0+0+mRowPrevious[1][j], 0+mRowPrevious[2][j])
			} else {
				mRowCurrent[0][j], trace[0][i][j] = tripleMaxTrace(scores[alpha[i-1]][beta[j-1]]+mRowPrevious[0][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[1][j-1], scores[alpha[i-1]][beta[j-1]]+mRowPrevious[2][j-1])
				mRowCurrent[1][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend+mRowCurrent[0][j-1], gapExtend+mRowCurrent[1][j-1], gapOpen+gapExtend+mRowCurrent[2][j-1])
				mRowCurrent[2][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend+mRowPrevious[0][j], gapOpen+gapExtend+mRowPrevious[1][j], gapExtend+mRowPrevious[2][j])
			}
		}
		if i < mColumn-1 {
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
		}
	}
	maxScore, route := affineTrace(mRowCurrent, mColumn, trace)
	return maxScore, route
}

func AffineGapChunk(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64, chunkSize int64) (int64, []Cigar) {
	var alphaSize, betaSize int64 = int64(len(alpha)), int64(len(beta))
	if alphaSize%chunkSize != 0 {
		log.Fatalf(fmt.Sprintf("Error: the first sequence, %s, has a length of %d, when it should be a multiple of %d\n", dna.BasesToString(alpha), alphaSize, chunkSize))
	}
	if betaSize%chunkSize != 0 {
		log.Fatalf(fmt.Sprintf("Error: the second sequence, %s, has a length of %d, when it should be a multiple of %d\n", dna.BasesToString(beta), betaSize, chunkSize))
	}
	alphaChunks := alphaSize / chunkSize
	betaChunks := betaSize / chunkSize

	mRowCurrent, mRowPrevious, mColumn, trace := initAffineScoringAndTrace(int(alphaChunks), int(betaChunks))

	var chunkScore int64
	for i := 0; i < mColumn; i++ {
		for j := range mRowCurrent[0] {
			if i == 0 && j == 0 {
				mRowCurrent[0][j] = 0
				mRowCurrent[1][j] = gapOpen
				mRowCurrent[2][j] = gapOpen
			} else if i == 0 {
				mRowCurrent[0][j] = veryNegNum
				mRowCurrent[1][j] = gapExtend*chunkSize + mRowCurrent[1][j-1]
				trace[1][i][j] = ColI
				mRowCurrent[2][j] = veryNegNum
			} else if j == 0 {
				mRowCurrent[0][j] = veryNegNum
				mRowCurrent[1][j] = veryNegNum
				mRowCurrent[2][j] = gapExtend*chunkSize + mRowPrevious[2][j]
				trace[2][i][j] = ColD
			} else {
				chunkScore = ungappedRegionScore(alpha, int64(i-1)*chunkSize, beta, int64(j-1)*chunkSize, chunkSize, scores)
				mRowCurrent[0][j], trace[0][i][j] = tripleMaxTrace(chunkScore+mRowPrevious[0][j-1], chunkScore+mRowPrevious[1][j-1], chunkScore+mRowPrevious[2][j-1])
				mRowCurrent[1][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend*chunkSize+mRowCurrent[0][j-1], gapExtend*chunkSize+mRowCurrent[1][j-1], gapOpen+gapExtend*chunkSize+mRowCurrent[2][j-1])
				mRowCurrent[2][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend*chunkSize+mRowPrevious[0][j], gapOpen+gapExtend*chunkSize+mRowPrevious[1][j], gapExtend*chunkSize+mRowPrevious[2][j])
			}
		}
		if i < mColumn-1 {
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
		}
	}
	maxScore, route := affineTrace(mRowCurrent, mColumn, trace)
	expandCigarRunLength(route, chunkSize)

	return maxScore, route
}

func multipleAffineGap(alpha []fasta.Fasta, beta []fasta.Fasta, scores [][]int64, gapOpen int64, gapExtend int64) (int64, []Cigar) {
	mRowCurrent, mRowPrevious, mColumn, trace := initAffineScoringAndTrace(len(alpha[0].Seq), len(beta[0].Seq))

	for i := 0; i < mColumn; i++ {
		for j := range mRowCurrent[0] {
			if i == 0 && j == 0 {
				mRowCurrent[0][j] = 0
				mRowCurrent[1][j] = gapOpen
				mRowCurrent[2][j] = gapOpen
			} else if i == 0 {
				mRowCurrent[0][j] = veryNegNum
				mRowCurrent[1][j] = gapExtend + mRowCurrent[1][j-1]
				trace[1][i][j] = ColI
				mRowCurrent[2][j] = veryNegNum
			} else if j == 0 {
				mRowCurrent[0][j] = veryNegNum
				mRowCurrent[1][j] = veryNegNum
				mRowCurrent[2][j] = gapExtend + mRowPrevious[2][j]
				trace[2][i][j] = ColD
			} else {
				mRowCurrent[0][j], trace[0][i][j] = tripleMaxTrace(scoreColumnMatch(alpha, beta, i-1, j-1, scores)+mRowPrevious[0][j-1], scoreColumnMatch(alpha, beta, i-1, j-1, scores)+mRowPrevious[1][j-1], scoreColumnMatch(alpha, beta, i-1, j-1, scores)+mRowPrevious[2][j-1])
				mRowCurrent[1][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend+mRowCurrent[0][j-1], gapExtend+mRowCurrent[1][j-1], gapOpen+gapExtend+mRowCurrent[2][j-1])
				mRowCurrent[2][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend+mRowPrevious[0][j], gapOpen+gapExtend+mRowPrevious[1][j], gapExtend+mRowPrevious[2][j])
			}
		}
		if i < mColumn-1 {
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
		}
	}
	maxScore, route := affineTrace(mRowCurrent, mColumn, trace)

	return maxScore, route
}

func multipleAffineGapChunk(alpha []fasta.Fasta, beta []fasta.Fasta, scores [][]int64, gapOpen int64, gapExtend int64, chunkSize int64) (int64, []Cigar) {
	var alphaSize, betaSize int64 = int64(len(alpha[0].Seq)), int64(len(beta[0].Seq))
	if alphaSize%chunkSize != 0 {
		log.Fatalf(fmt.Sprintf("Error: the first subalignment has a length of %d, when it should be a multiple of %d\n", alphaSize, chunkSize))
	}
	if betaSize%chunkSize != 0 {
		log.Fatalf(fmt.Sprintf("Error: the second subalignment has a length of %d, when it should be a multiple of %d\n", betaSize, chunkSize))
	}
	alphaChunks := alphaSize / chunkSize
	betaChunks := betaSize / chunkSize

	mRowCurrent, mRowPrevious, mColumn, trace := initAffineScoringAndTrace(int(alphaChunks), int(betaChunks))

	var chunkScore int64
	for i := 0; i < mColumn; i++ {
		for j := range mRowCurrent[0] {
			if i == 0 && j == 0 {
				mRowCurrent[0][j] = 0
				mRowCurrent[1][j] = gapOpen
				mRowCurrent[2][j] = gapOpen
			} else if i == 0 {
				mRowCurrent[0][j] = veryNegNum
				mRowCurrent[1][j] = gapExtend*chunkSize + mRowCurrent[1][j-1]
				trace[1][i][j] = ColI
				mRowCurrent[2][j] = veryNegNum
			} else if j == 0 {
				mRowCurrent[0][j] = veryNegNum
				mRowCurrent[1][j] = veryNegNum
				mRowCurrent[2][j] = gapExtend*chunkSize + mRowPrevious[2][j]
				trace[2][i][j] = ColD
			} else {
				chunkScore = ungappedRegionColumnScore(alpha, (i-1)*int(chunkSize), beta, (j-1)*int(chunkSize), int(chunkSize), scores)
				mRowCurrent[0][j], trace[0][i][j] = tripleMaxTrace(chunkScore+mRowPrevious[0][j-1], chunkScore+mRowPrevious[1][j-1], chunkScore+mRowPrevious[2][j-1])
				mRowCurrent[1][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend*chunkSize+mRowCurrent[0][j-1], gapExtend*chunkSize+mRowCurrent[1][j-1], gapOpen+gapExtend*chunkSize+mRowCurrent[2][j-1])
				mRowCurrent[2][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend*chunkSize+mRowPrevious[0][j], gapOpen+gapExtend*chunkSize+mRowPrevious[1][j], gapExtend*chunkSize+mRowPrevious[2][j])
			}
		}
		if i < mColumn-1 {
			mRowPrevious, mRowCurrent = mRowCurrent, mRowPrevious
		}
	}
	maxScore, route := affineTrace(mRowCurrent, mColumn, trace)
	expandCigarRunLength(route, chunkSize)

	return maxScore, route
}

func scoreAffineAln(alpha fasta.Fasta, beta fasta.Fasta, scores [][]int64, gapOpen int64, gapExtend int64) (int64, error) {
	if len(alpha.Seq) != len(beta.Seq) {
		return 0, fmt.Errorf("Error: alignment being scored has sequences of unequal length: %d, %d\n", len(alpha.Seq), len(beta.Seq))
	}
	var score int64 = 0
	alphaInGap, betaInGap := false, false
	for i := range alpha.Seq {
		if alpha.Seq[i] != dna.Gap && beta.Seq[i] != dna.Gap {
			score += scores[alpha.Seq[i]][beta.Seq[i]]
		}
		if alpha.Seq[i] == dna.Gap {
			if !alphaInGap {
				score += gapOpen
			}
			score += gapExtend
			alphaInGap = true
		} else {
			alphaInGap = false
		}
		if beta.Seq[i] == dna.Gap {
			if !betaInGap {
				score += gapOpen
			}
			score += gapExtend
			betaInGap = true
		} else {
			betaInGap = false
		}
	}
	return score, nil
}
