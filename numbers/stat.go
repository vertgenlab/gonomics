package numbers

import (
	"math"
)

func AverageFloat64(f []float64) float64 {
	answer := 0.0
	for i := 0; i < len(f); i++ {
		answer = answer + f[i]
	}
	return answer / float64(len(f))
}

func VarianceFloat64(f []float64) float64 {
	answer := 0.0
	average := AverageFloat64(f)
	for i := 0; i < len(f); i++ {
		answer = answer + math.Pow(f[i]-average, 2)
	}
	return answer / float64(len(f)-1)
}

func StandardDeviationFloat64(f []float64) float64 {
	return math.Sqrt(VarianceFloat64(f))
}
