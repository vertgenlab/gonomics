package popgen

import (
	"github.com/vertgenlab/gonomics/numbers"
	"math"
)

//LogIntegrateStationarityCache is a specialized numerical integration function for logspace integrals over the stationarity distribution, used for MCMC.
func LogIntegrateStationarityCache(alpha float64, n int, k int, stepSize int, alleleFrequencyCache []float64, nkpCache [][][]float64) float64 {
	var bins int = len(nkpCache[n][k])-1
	var deltaX float64 = (float64(stepSize) * (RightBound - LeftBound)) / float64(bins)
	var logDeltaX float64 = math.Log(deltaX)
	var currLeftIndex int = 0
	var currRightIndex int = stepSize
	var sumHeights, rightEval float64

	//first tiem, sets answer as the area of the first rectangle
	var nextLeftEval float64 = f(alpha, n, k, currRightIndex, alleleFrequencyCache, nkpCache)
	var firstLeft float64 = f(alpha, n, k, currLeftIndex, alleleFrequencyCache, nkpCache)
	sumHeights = numbers.MidpointLog(firstLeft, nextLeftEval)

	for i := stepSize; i < bins; i += stepSize {
		currLeftIndex += stepSize
		currRightIndex += stepSize
		rightEval = f(alpha, n, k, currRightIndex, alleleFrequencyCache, nkpCache)
		rightEval = f(alpha, n, k, currRightIndex, alleleFrequencyCache, nkpCache)
		sumHeights = numbers.AddLog(sumHeights, numbers.MidpointLog(nextLeftEval, rightEval))
		nextLeftEval = rightEval
	}
	return numbers.MultiplyLog(sumHeights, logDeltaX)
}

//helper function of LogIntegrateStationarityCache, returns the function evaluation for the unadjusted stationarity distribution.
func f(alpha float64, n int, k int, index int, alleleFrequencyCache []float64, nkpCache [][][]float64) float64 {
	return numbers.MultiplyLog(math.Log(AFSStationarity(alleleFrequencyCache[index], alpha)), nkpCache[n][k][index])
}

//helper function for LogintegrateStationarityCache,  returns the function evaluation for the stationarity distribution adjusted for rare allele detection bias.
func fCorrected(alpha float64, n int, k int, index int, alleleFrequencyCache []float64, nkpCache [][][]float64) float64 {
	return numbers.MultiplyLog(f(alpha, n, k, index, alleleFrequencyCache, nkpCache), DetectionProbability(index, n, nkpCache))
}
