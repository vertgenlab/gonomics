package align

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

// the data structure is a 3d slice where the first index is 0,1,2 and represents
// the match, gap in x (first seq), and gap in y (second seq).
func initAffineScoringAndTrace(firstSeqLen int, secondSeqLen int) ([][][]int64, [][][]ColType) {
	m := make([][][]int64, 3)
	trace := make([][][]ColType, 3)
	for k := range m {
		m[k] = make([][]int64, firstSeqLen+1)
		trace[k] = make([][]ColType, firstSeqLen+1)
		for i := range m[0] {
			m[k][i] = make([]int64, secondSeqLen+1)
			trace[k][i] = make([]ColType, secondSeqLen+1)
		}
	}
	return m, trace
}

func affineTrace(m [][][]int64, trace [][][]ColType) (int64, []Cigar) {
	route := make([]Cigar, 1)
	lastI := len(m[0]) - 1
	lastJ := len(m[0][0]) - 1
	maxScore, k := tripleMaxTrace(m[0][lastI][lastJ], m[1][lastI][lastJ], m[2][lastI][lastJ])
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

func AffineGap(alpha []dna.Base, beta []dna.Base, scores [][]int64, gapOpen int64, gapExtend int64) (int64, []Cigar) {
	m, trace := initAffineScoringAndTrace(len(alpha), len(beta))
	for i := range m[0] {
		for j := range m[0][0] {
			if i == 0 && j == 0 {
				m[0][i][j] = 0
				m[1][i][j] = gapOpen
				m[2][i][j] = gapOpen
			} else if i == 0 {
				m[0][i][j] = veryNegNum
				m[1][i][j] = gapExtend + m[1][i][j-1]
				trace[1][i][j] = ColI /*new*/
				m[2][i][j] = veryNegNum
			} else if j == 0 {
				m[0][i][j] = veryNegNum
				m[1][i][j] = veryNegNum
				m[2][i][j] = gapExtend + m[2][i-1][j]
				trace[2][i][j] = ColD /*new*/
			} else {
				m[0][i][j], trace[0][i][j] = tripleMaxTrace(scores[alpha[i-1]][beta[j-1]]+m[0][i-1][j-1], scores[alpha[i-1]][beta[j-1]]+m[1][i-1][j-1], scores[alpha[i-1]][beta[j-1]]+m[2][i-1][j-1])
				m[1][i][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend+m[0][i][j-1], gapExtend+m[1][i][j-1], gapOpen+gapExtend+m[2][i][j-1])
				m[2][i][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend+m[0][i-1][j], gapOpen+gapExtend+m[1][i-1][j], gapExtend+m[2][i-1][j])
			}
		}
	}
	maxScore, route := affineTrace(m, trace)

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

	m, trace := initAffineScoringAndTrace(int(alphaChunks), int(betaChunks))

	var chunkScore int64
	for i := range m[0] {
		for j := range m[0][0] {
			if i == 0 && j == 0 {
				m[0][i][j] = 0
				m[1][i][j] = gapOpen
				m[2][i][j] = gapOpen
			} else if i == 0 {
				m[0][i][j] = veryNegNum
				m[1][i][j] = gapExtend*chunkSize + m[1][i][j-1]
				trace[1][i][j] = ColI
				m[2][i][j] = veryNegNum
			} else if j == 0 {
				m[0][i][j] = veryNegNum
				m[1][i][j] = veryNegNum
				m[2][i][j] = gapExtend*chunkSize + m[2][i-1][j]
				trace[2][i][j] = ColD
			} else {
				chunkScore = ungappedRegionScore(alpha, int64(i-1)*chunkSize, beta, int64(j-1)*chunkSize, chunkSize, scores)
				m[0][i][j], trace[0][i][j] = tripleMaxTrace(chunkScore+m[0][i-1][j-1], chunkScore+m[1][i-1][j-1], chunkScore+m[2][i-1][j-1])
				m[1][i][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend*chunkSize+m[0][i][j-1], gapExtend*chunkSize+m[1][i][j-1], gapOpen+gapExtend*chunkSize+m[2][i][j-1])
				m[2][i][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend*chunkSize+m[0][i-1][j], gapOpen+gapExtend*chunkSize+m[1][i-1][j], gapExtend*chunkSize+m[2][i-1][j])
			}
		}
	}

	maxScore, route := affineTrace(m, trace)
	expandCigarRunLength(route, chunkSize)

	return maxScore, route
}

func multipleAffineGap(alpha []fasta.Fasta, beta []fasta.Fasta, scores [][]int64, gapOpen int64, gapExtend int64) (int64, []Cigar) {
	m, trace := initAffineScoringAndTrace(len(alpha[0].Seq), len(beta[0].Seq))

	for i := range m[0] {
		for j := range m[0][0] {
			if i == 0 && j == 0 {
				m[0][i][j] = 0
				m[1][i][j] = gapOpen
				m[2][i][j] = gapOpen
			} else if i == 0 {
				m[0][i][j] = veryNegNum
				m[1][i][j] = gapExtend + m[1][i][j-1]
				trace[1][i][j] = ColI
				m[2][i][j] = veryNegNum
			} else if j == 0 {
				m[0][i][j] = veryNegNum
				m[1][i][j] = veryNegNum
				m[2][i][j] = gapExtend + m[2][i-1][j]
				trace[2][i][j] = ColD
			} else {
				m[0][i][j], trace[0][i][j] = tripleMaxTrace(scoreColumnMatch(alpha, beta, i-1, j-1, scores)+m[0][i-1][j-1], scoreColumnMatch(alpha, beta, i-1, j-1, scores)+m[1][i-1][j-1], scoreColumnMatch(alpha, beta, i-1, j-1, scores)+m[2][i-1][j-1])
				m[1][i][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend+m[0][i][j-1], gapExtend+m[1][i][j-1], gapOpen+gapExtend+m[2][i][j-1])
				m[2][i][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend+m[0][i-1][j], gapOpen+gapExtend+m[1][i-1][j], gapExtend+m[2][i-1][j])
			}
		}
	}

	maxScore, route := affineTrace(m, trace)
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

	m, trace := initAffineScoringAndTrace(int(alphaChunks), int(betaChunks))

	var chunkScore int64
	for i := range m[0] {
		for j := range m[0][0] {
			if i == 0 && j == 0 {
				m[0][i][j] = 0
				m[1][i][j] = gapOpen
				m[2][i][j] = gapOpen
			} else if i == 0 {
				m[0][i][j] = veryNegNum
				m[1][i][j] = gapExtend*chunkSize + m[1][i][j-1]
				trace[1][i][j] = ColI
				m[2][i][j] = veryNegNum
			} else if j == 0 {
				m[0][i][j] = veryNegNum
				m[1][i][j] = veryNegNum
				m[2][i][j] = gapExtend*chunkSize + m[2][i-1][j]
				trace[2][i][j] = ColD
			} else {
				chunkScore = ungappedRegionColumnScore(alpha, (i-1)*int(chunkSize), beta, (j-1)*int(chunkSize), int(chunkSize), scores)
				m[0][i][j], trace[0][i][j] = tripleMaxTrace(chunkScore+m[0][i-1][j-1], chunkScore+m[1][i-1][j-1], chunkScore+m[2][i-1][j-1])
				m[1][i][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend*chunkSize+m[0][i][j-1], gapExtend*chunkSize+m[1][i][j-1], gapOpen+gapExtend*chunkSize+m[2][i][j-1])
				m[2][i][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend*chunkSize+m[0][i-1][j], gapOpen+gapExtend*chunkSize+m[1][i-1][j], gapExtend*chunkSize+m[2][i-1][j])
			}
		}
	}

	maxScore, route := affineTrace(m, trace)
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
