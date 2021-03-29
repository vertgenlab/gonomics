package numbers

import (
	"math"
)

//AverageFloat64 returns the mean value of type float64 from a slice of type float64.
func AverageFloat64(f []float64) float64 {
	answer := 0.0
	for i := range f {
		answer = answer + f[i]
	}
	return answer / float64(len(f))
}

//VarianceFloat64 returns the variance (type float64) from a slice of type float64.
func VarianceFloat64(f []float64) float64 {
	answer := 0.0
	average := AverageFloat64(f)
	for i := range f {
		answer = answer + math.Pow(f[i]-average, 2)
	}
	return answer / float64(len(f)-1)
}

//StandardDeviationFloat64 returns the standard deviation (type float64) from a slice of type float64.
func StandardDeviationFloat64(f []float64) float64 {
	return math.Sqrt(VarianceFloat64(f))
}
