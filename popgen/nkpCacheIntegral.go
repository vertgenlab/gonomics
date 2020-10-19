package popgen

import (
	"github.com/vertgenlab/gonomics/numbers"
	"math"
)

func LogIntegrateStationarityCache(alpha float64, n int, k int, stepSize int, alleleFrequencyCache []float64, nkpCache [][][]float64) float64 {
	var bins int = len(nkpCache[n][k])-1
	var deltaX float64 = (RightBound - LeftBound) / float64(bins)
	var logDeltaX float64 = math.Log(deltaX)
	var currLeftIndex int = 0
	var currRightIndex int = stepSize
	var answer, rightEval float64

	//first tiem, sets answer as the area of the first rectangle
	var nextLeftEval float64 = math.Log(AFSStationarity(alleleFrequencyCache[currRightIndex], nkpCache[n][k][currRightIndex]))
	answer = f(alpha, currLeftIndex, alleleFrequencyCache, nkpCache[n][k][currLeftIndex])
	firstLeft := math.Log(AFSStationarity(alleleFrequencyCache[currLeftIndex], nkpCache[n][k][currLeftIndex]))
	answer = numbers.MultiplyLog(numbers.MidpointLog(firstLeft, nextLeftEval), logDeltaX)

	for i := stepSize; i < bins; i += stepSize {
		currLeftIndex += stepSize
		currRightIndex += stepSize
		rightEval = f(alpha, currRightIndex, alleleFrequencyCache, nkpCache[n][k][currRightIndex])
		answer = numbers.AddLog(answer, numbers.MultiplyLog(numbers.MidpointLog(nextLeftEval, rightEval), logDeltaX))
		nextLeftEval = rightEval
	}
	return answer
}

func f(alpha float64, index int, alleleFrequencyCache[]float64, binomial float64) float64 {
	return numbers.MultiplyLog(math.Log(AFSStationarity(alleleFrequencyCache[index], alpha)), binomial)
}
