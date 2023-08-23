// package fit contains tools for fitting statistical distributions to data.

package fit

import (
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"log"
	"math"
	"math/rand"
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

func zeroInflatedNegativeBinomialLogLikelihood(data []int, R float64, P float64, Z float64) float64 {
	likelihood := 0.0
	var density float64
	for i := range data {
		if i == 0 {
			likelihood += float64(data[0]) * math.Log(Z+(1-Z)*math.Pow(P, R))
			//fmt.Printf("AtZero: R: %v. P: %v. Z: %v. Likelihood: %v.\n", R, P, Z, likelihood)
		} else {
			density, _ = numbers.NegativeBinomialDist(i, R, P, true)
			likelihood += float64(data[i]) * logspace.Multiply(math.Log(1-Z), density)
			//fmt.Printf("AtGreaterThanZero: R: %v. P: %v. Z: %v. Likelihood: %v.\n", R, P, Z, likelihood)
		}
	}
	return likelihood
}

func ZeroInflatedNegativeBinomial(data []int, learningRate float64, delta float64, epsilon float64) (float64, float64, float64) {
	var R, P, Z float64 = 1.0, 0.5, 0.01 //hardcoded initialization
	var prevR, prevP, prevZ = 1.0, 0.5, 0.01
	//var R, P, Z float64 = 2.0, 0.3, 0.2 //hardcoded initialization
	var gradient = make([]float64, 3) // gradient has three partial derivatives, one per model parameter
	var lossPlus, lossMinus float64

	var propDiff float64 = 1 // proportion difference between current loss and previous loss
	var currLoss float64
	var prevLoss float64 = zeroInflatedNegativeBinomialLogLikelihood(data, R, P, Z)

	for propDiff > epsilon {
		R += delta
		lossPlus = zeroInflatedNegativeBinomialLogLikelihood(data, R, P, Z)
		R -= 2 * delta
		lossMinus = zeroInflatedNegativeBinomialLogLikelihood(data, R, P, Z)
		R += delta
		gradient[0] = (lossPlus - lossMinus) / (2 * delta)

		P += delta
		lossPlus = zeroInflatedNegativeBinomialLogLikelihood(data, R, P, Z)
		P -= 2 * delta
		lossMinus = zeroInflatedNegativeBinomialLogLikelihood(data, R, P, Z)
		P += delta
		gradient[1] = (lossPlus - lossMinus) / (2 * delta)

		Z += delta
		lossPlus = zeroInflatedNegativeBinomialLogLikelihood(data, R, P, Z)
		Z -= 2 * delta
		lossMinus = zeroInflatedNegativeBinomialLogLikelihood(data, R, P, Z)
		Z += delta
		gradient[2] = (lossPlus - lossMinus) / (2 * delta)

		//fmt.Printf("GradR:\t%v\tGradP:\t%v\tGradZ:\t%v\n", gradient[0], gradient[1], gradient[2])

		prevR, prevP, prevZ = R, P, Z

		if learningRate*gradient[0] > 0.1 {
			R += 0.1
			if R < 0 {
				R -= 0.1
			}
		} else if learningRate*gradient[0] < -0.1 {
			R += -0.1
			if R < 0 {
				R -= -0.1
			}
		} else {
			R += learningRate * gradient[0]
			if R < 0 {
				R -= learningRate * gradient[0]
			}
		}

		if learningRate*gradient[1] > 0.1 {
			P += 0.1
			if P >= 1 || P <= 0 {
				P -= 0.1
			}
		} else if learningRate*gradient[1] < -0.1 {
			P += -0.1
			if P > 1 || P < 0 {
				P -= -0.1
			}
		} else {
			P += learningRate * gradient[1]
			if P > 1 || P < 0 {
				P -= learningRate * gradient[1]
			}
		}

		if learningRate*gradient[2] > 0.1 {
			Z += 0.1
			if Z < 0 || Z > 1 {
				Z -= 0.1
			}
		} else if learningRate*gradient[2] < -0.1 {
			Z -= 0.1
			if Z < 0 || Z > 1 {
				Z -= -0.1
			}
		} else {
			Z += learningRate * gradient[2]
			if Z < 0 || Z > 1 {
				Z -= learningRate * gradient[2]
			}
		}

		currLoss = zeroInflatedNegativeBinomialLogLikelihood(data, R, P, Z)
		//fmt.Printf("LossFunction:\t%v\tR:\t%v\tP:\t%v\tZ:\t%v\n", currLoss, R, P, Z)

		if currLoss < prevLoss {
			return prevR, prevP, prevZ
		}
		propDiff = math.Abs(math.Abs(currLoss-prevLoss) / prevLoss)
		prevLoss = currLoss
	}

	return R, P, Z
}

func simulateZeroInflatedNegativeBinomial(r float64, p float64, zeroExcessParameter float64, samples int) []int {
	var data = make([]int, 10)
	var tmpData []int
	var currVariate int
	if r < 1 {
		log.Fatalf("Error: Cannot simulate zero-inflated negative binomial with r < 1. Found: %v\n", r)
	}
	if p < 0 || p > 1 {
		log.Fatalf("Error: to simulate zero-inflated negative binomial, p must be a value between zero and one. Found: %v.\n", p)
	}
	if zeroExcessParameter > 1 || zeroExcessParameter < 0 {
		log.Fatalf("Error: the zeroExcessParameter must be a value between zero and one. Found: %v\n", zeroExcessParameter)
	}
	for i := 0; i < samples; i++ {
		currVariate = simulateZeroInflatedNegativeBinomialValue(r, p, zeroExcessParameter)
		if currVariate+1 > len(data) {
			tmpData = make([]int, currVariate+1)
			copy(tmpData, data)
			data = tmpData
		}
		data[currVariate]++
	}
	return data
}

func simulateZeroInflatedNegativeBinomialValue(r float64, p float64, zeroExcessParameter float64) int {
	if rand.Float64() < zeroExcessParameter {
		return 0
	}
	return randNegativeBinomial(r, p)
}

func randNegativeBinomial(r float64, p float64) int {
	s := 0           //successes
	f := 0           //failures
	for s < int(r) { // for fewer than r successes
		if rand.Float64() < p {
			s++
		} else {
			f++
		}
	}
	return f
}
