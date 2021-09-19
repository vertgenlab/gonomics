package align

import (
	"github.com/vertgenlab/gonomics/dna"
)

func ungappedRegionScore(alpha []dna.Base, alphaStart int, beta []dna.Base, betaStart int, length int, scores [][]int) int {
	var answer int = 0
	for i, j := alphaStart, betaStart; i < alphaStart+length; i, j = i+1, j+1 {
		answer += scores[alpha[i]][beta[j]]
	}
	return answer
}
