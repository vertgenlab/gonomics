package simulate

import (
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math"
)

func getStation(n int, k int, alpha float64) func(float64) float64 {
	var binomCoeff float64 = numbers.BinomCoefficientLog(n, k)
	return func(p float64) float64 {
		expression := numbers.BinomialExpressionLog(n-2, k-1, p)
		logPart := math.Log((1.0 - math.Exp(-1.0*alpha*(1.0-p))) * 2.0 / (1.0 - math.Exp(-1.0*alpha)))
		return numbers.MultiplyLog(binomCoeff, numbers.MultiplyLog(expression, logPart))
	}
}

func binomInIntegral(n int, k int, alpha float64, start float64, end float64, accuracy float64) float64 {
	f := getStation(n, k, alpha)
	integralResult := numbers.AdaptiveSimpsonsLog(f, start, end, accuracy, 1000)
	return integralResult
}

func binomInIntegralParts(n int, k int, alpha float64, accuracy float64) float64 {
	var switchPoint float64 = float64(k) / float64(n)
	if switchPoint <= 0 || switchPoint >= 1 {
		log.Fatal("Error: confused about switch point\n")
	}
	f := getStation(n, k, alpha)
	partOne := numbers.AdaptiveSimpsonsLog(f, 0, switchPoint, accuracy, 1000)
	partTwo := numbers.AdaptiveSimpsonsLog(f, switchPoint, 1, accuracy, 1000)
	return numbers.AddLog(partOne, partTwo)
}

func AFSLikelihoodNew(count []int, totalAlleles int, alpha float64, accuracy float64) float64 {
	var answer float64 = 0.0
	for j := 0; j < len(count); j++ {
		answer = numbers.MultiplyLog(answer, AlleleFrequencyProbabilityNew(totalAlleles, count[j], alpha, accuracy))
	}
	return answer
}

func AlleleFrequencyProbabilityNew(n int, i int, alpha float64, accuracy float64) float64 {
	var denom float64 = math.Inf(-1)
	var numer float64

	for j := 1; j < n; j++ {
		denom = numbers.AddLog(denom, binomInIntegralParts(n, j, alpha, accuracy))
	}
	numer = binomInIntegralParts(n, i, alpha, accuracy) // original answer log, integral log

	return numbers.DivideLog(numer, denom)
}
