package numbers

import (
	"log"
	"math"
)

// AverageFloat64 returns the mean value of type float64 from a slice of type float64.
func AverageFloat64(f []float64) float64 {
	answer := 0.0
	for i := range f {
		answer = answer + f[i]
	}
	return answer / float64(len(f))
}

// VarianceFloat64 returns the variance (type float64) from a slice of type float64.
func VarianceFloat64(f []float64) float64 {
	answer := 0.0
	average := AverageFloat64(f)
	for i := range f {
		answer = answer + math.Pow(f[i]-average, 2)
	}
	return answer / float64(len(f)-1)
}

// StandardDeviationFloat64 returns the standard deviation (type float64) from a slice of type float64.
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

// BenjaminiHochberg takes in a []float64 that contains sorted pValues from smallest to largest. It returns
// Benjamini-Hochberg adjusted pValues according to the formula pAdj = p (n/k) where p is the unadjusted pValue,
// k is the pValue's rank and n is the total number of p values in the slice. If pValues are not sorted, the program will log Fatal
func BenjaminiHochberg(pVals [][]float64) [][]float64 {
	var adjP [][]float64
	var newP float64
	n := len(pVals)
	prevP := pVals[len(pVals)-1][1]
	for i := len(pVals) - 1; i >= 0; i-- {
		if i > 0 {
			if pVals[i][1] < pVals[i-1][1] {
				log.Fatalf("input pValues must be sorted from smallest to largest")
			}
		}
		newP = pVals[i][1] * (float64(n) / float64(i+1))
		if newP >= prevP {
			adjP = append([][]float64{{pVals[i][0], prevP}}, adjP...)
		} else {
			adjP = append([][]float64{{pVals[i][0], newP}}, adjP...)
			prevP = newP
		}
	}
	return adjP
}
