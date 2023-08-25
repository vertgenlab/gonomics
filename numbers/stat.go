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

func AverageInt(d []int) float64 {
	answer := 0
	for _, i := range d {
		answer = answer + i
	}
	return float64(answer) / float64(len(d))
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

/*
func TwoSideTTest(a []float64, b []float64) float64 {
	aAvg := AverageFloat64(a)
	bAvg := AverageFloat64(b)
	aSD := StandardDeviationFloat64(a)
	bSD := StandardDeviationFloat64(b)

	numerator := aAvg - bAvg
	d := (aSD*aSD)/float64(len(a)) + (bSD*bSD)/float64(len(b))

	t := numerator / math.Sqrt(d)
	fmt.Println("t: ", t)
	dof := len(a) + len(b) - 2
	return cumulativeDistributionFunction(t, dof)
}

func cumulativeDistributionFunction(t float64, dof int) float64 {
	term1 := 0.5 + t*math.Gamma(float64(dof+1)/float64(2))
	b := float64(dof+1) / float64(2)
	c := float64(3) / float64(2)
	z := -math.Pow(t, 2) / float64(dof)
	term2Num := HypergeometricFunction(0.5, b, c, z)
	term2Denom := math.Sqrt(math.Pi*float64(dof)) * math.Gamma(float64(dof)/2)
	return term1 * (term2Num / term2Denom)
}*/
