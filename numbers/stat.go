package numbers

import (
	"log"
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

// Pearson calculates the Pearson Correlation Coefficient between two slices of float64.
func Pearson(a []float64, b []float64) float64 {
	if len(a) != len(b) {
		log.Fatalf("Error in Pearson. Input float slices must be of the same length.")
	}
	if len(a) == 0 {
		log.Fatalf("Cannot compute the Pearson Correlation Coefficient for empty vectors.")
	}
	aAve := AverageFloat64(a)
	bAve := AverageFloat64(b)
	var num float64 //numerator
	var sumA, sumB float64
	for i := range a {
		num += (a[i] - aAve) * (b[i] - bAve)
		sumA += (a[i] - aAve) * (a[i] - aAve)
		sumB += (b[i] - bAve) * (b[i] - bAve)
	}
	return num / (math.Sqrt(sumA) * math.Sqrt(sumB))
}
