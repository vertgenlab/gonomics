package align

import (
	"github.com/craiglowe/gonomics/dna"
)

func ungappedRegionScore(alpha []dna.Base, alphaStart int64, beta []dna.Base, betaStart int64, length int64, scores [][]int64) int64 {
	var answer int64 = 0
	for i, j := alphaStart, betaStart; i < alphaStart+length; i, j = i+1, j+1 {
		answer += scores[alpha[i]][beta[j]]
	}
	return answer
}
