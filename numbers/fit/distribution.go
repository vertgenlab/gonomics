// package fit contains tools for fitting statistical distributions to data.

package fit

import (
	"github.com/vertgenlab/gonomics/numbers"
)

// Poisson fits an input slice of data points []float64 to a Poisson distribution, defined
// by the single parameter lambda, which is returned as a float64.
func Poisson(data []float64) float64 {
	return numbers.AverageFloat64(data)
}

// PoissonHistogram fits a histogram to a Poisson distribution, defined by the single
// parameter lambda, returned as a float64.
// the histogram is formatted by index:value pairs which correspond to the observed number (index)
// and the number of times that observation was seen (value).
func PoissonHistogram(histogram []int) float64 {
	totalNumOfBins := 0
	sumOfAllBinValues := 0
	for bin := 0; bin < len(histogram); bin++ {
		totalNumOfBins += histogram[bin]
		sumOfAllBinValues += bin * histogram[bin]
	}
	return float64(sumOfAllBinValues) / float64(totalNumOfBins)
}

// NegativeBinomial fits an input slice of data points []float64 to a negative binomial distribution,
// defined by two float64 parameters, which are returned.
// the first return is 'r', also called "shape".
// the second return is 'p', or the success probability.
// The third bool returns true when a distribution could not be fit to the input mean and variance.
func NegativeBinomial(data []float64) (float64, float64, bool) {
	average := numbers.AverageFloat64(data)
	variance := numbers.VarianceFloat64(data)
	return NegativeBinomialFromSumStats(average, variance)
}

// NegativeBinomialFromSumStats fits a negative binomial distribution to a dataset, represented
// by the summary statistics 'mean' and 'variance'.
// the first return is 'r', also called "shape".
// the second return is 'p', or the success probability.
// The third bool returns true when a distribution could not be fit to the input mean and variance.
func NegativeBinomialFromSumStats(mean float64, variance float64) (float64, float64, bool) {
	if variance <= 0 || mean <= 0 {
		return -1, -1, true
	}
	p := mean / variance
	if p <= 0 || p >= 1 {
		return -1, -1, true
	}

	r := (mean * p) / (1 - p)
	if r < 0 {
		return -1, -1, true
	}

	return r, p, false
}

func NegativeBinomialFromCountSlice(data []int) (float64, float64, bool) {
	records := make([]float64, 0)
	for i := range data {
		for j := 0; j < data[i]; j++ {
			records = append(records, float64(i))
		}
	}

	return NegativeBinomial(records)
}
