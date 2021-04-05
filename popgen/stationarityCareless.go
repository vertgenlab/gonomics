package popgen

import (
	"github.com/vertgenlab/gonomics/numbers"
	"math"
)

//TODO: (Riley) This file has a few old versions of functions that I'm keeping around while debugging. Will eventually delete.

//FIntegralComponentCareless is an outdated version of FIntegralComponent. As the new version does a "careful" calculation, this one is "careless". Used for plotting and presentation figures.
func FIntegralComponentCareless(n int, k int, alpha float64) func(float64) float64 {
	var binomCoeff float64 = numbers.BinomCoefficientLog(n, k)
	return func(p float64) float64 {
		carefulExpression := numbers.MultiplyLog(numbers.LogPow(p, float64(k-1)), numbers.LogPow(1.0-p, float64(n-k-1)))
		var ans float64 = binomCoeff + carefulExpression + math.Log((1.0-math.Exp(-1.0*alpha*(1.0-p)))*2.0/(1.0-math.Exp(-1.0*alpha)))
		return math.Exp(ans)
	}
}

//AfsSampleDensityCareless is an outdated version of AfsSampleDensity. As the new version does a "careful" calculation, this one is "careless". Used for plotting and presentation figures.
func AfsSampleDensityCareless(n int, k int, alpha float64) float64 {
	var switchPoint float64 = float64(k) / float64(n)
	f := FIntegralComponentCareless(n, k, alpha)
	integral := numbers.AddLog(numbers.AdaptiveSimpsonsLog(f, 0.0, switchPoint, 1e-8, 100000), numbers.AdaptiveSimpsonsLog(f, switchPoint, 1.0, 1e-8, 100))
	return integral
}
