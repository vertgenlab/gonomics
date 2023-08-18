// package fit contains tools for fitting statistical distributions to data.

package fit

import (
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

// Poisson fits an input slice of data points []float64 to a Poisson distribution, defined
// by the single parameter lambda, which is returned as a float64.
func Poisson(data []float64) float64 {
	return numbers.AverageFloat64(data)
}

// NegativeBinomial fits an input slice of data points []float64 to a negative binomial distribution,
// defined by two float64 parameters, which are returned.
// the first return is 'r', also called "shape".
// the second return is 'p', or the success probability.
func NegativeBinomial(data []float64) (float64, float64) {
	average := numbers.AverageFloat64(data)
	variance := numbers.VarianceFloat64(data)
	if variance <= 0 {
		log.Fatalf("Error: cannot fit negative binomial to data, as the variance must be positive.")
	}
	p := average / variance
	if p <= 0 {
		log.Fatalf("Error: cannot fit negative binomial to data, as p (mean / variance) must be positive. ")
	}
	r := (average * p) / (1 - p)
	return r, p
}
