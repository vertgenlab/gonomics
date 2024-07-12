package numbers

import (
	"log"
	"math"
	"math/rand"
)

// SampleInverseNormal returns a simulated value from a normal distribution.
func SampleInverseNormal(mu float64, sigma float64, seed *rand.Rand) float64 {
	return seed.NormFloat64()*sigma + mu
}

// InitializeFastRejectionSampler takes in the parameters of a rejection sampler and returns the binHeights and sumHeights variables.
func InitializeFastRejectionSampler(xLeft float64, xRight float64, f func(float64) float64, bins int) ([]float64, float64) {
	if xLeft >= xRight {
		log.Fatalf("Error in FastRejectionSample: xRight must be greater than xLeft.")
	}
	var binHeights []float64 = make([]float64, bins)
	var sumHeights float64
	var support float64 = xRight - xLeft
	var stepSize float64 = support / float64(bins)
	var currLeft float64 = xLeft
	var currRight float64 = xLeft + stepSize
	var fCurrLeft, fCurrRight float64
	var firstTime bool = true

	for i := 0; i < bins; i++ {
		if firstTime {
			firstTime = false
			fCurrLeft = f(currLeft)
			fCurrRight = f(currRight)
			binHeights[i] = Max(fCurrLeft, fCurrRight)
		} else {
			fCurrLeft = fCurrRight
			currRight += stepSize
			fCurrRight = f(currRight)
			binHeights[i] = Max(fCurrLeft, fCurrRight)
		}
		sumHeights += binHeights[i]
	}
	return binHeights, sumHeights
}

// FastRejectionSampler returns simulated values from an a func(float64) float64 between a left and right value using an optimized rejection sampler
// that divides the function support into discrete bins with optimized sampling heights.
// maxSampleDepth triggers the log.Fatalf in the RejectionSample func, and samples is the number of values to be returned.
func FastRejectionSampler(xLeft float64, xRight float64, f func(float64) float64, bins int, maxSampleDepth int, samples int, seed *rand.Rand) []float64 {
	var answer []float64 = make([]float64, samples)
	var stepSize float64 = (xRight - xLeft) / float64(bins)
	binHeights, sumHeights := InitializeFastRejectionSampler(xLeft, xRight, f, bins)
	for j := 0; j < samples; j++ {
		answer[j] = RejectionSampleChooseBin(xLeft, xRight, stepSize, f, maxSampleDepth, sumHeights, binHeights, seed)
	}
	return answer
}

// RejectionSampleChooseBin is a helper function of FAstRejectionSampler.
func RejectionSampleChooseBin(xLeft float64, xRight float64, stepSize float64, f func(float64) float64, maxIteration int, sumHeights float64, binHeights []float64, seed *rand.Rand) float64 {
	var x, y float64
	var currBin int
	var currLeft, currRight float64
	for i := 0; i < maxIteration; i++ {
		currBin = chooseBin(sumHeights, binHeights, seed)
		currLeft = xLeft + float64(currBin)*stepSize
		currRight = currLeft + stepSize
		x = RandFloat64InRange(currLeft, currRight, seed)
		y = f(x)
		if RandFloat64InRange(0.0, binHeights[currBin], seed) < y {
			return x
		}
	}
	log.Fatalf("Exceeded max iteration in RejectionSampleChooseBin.")
	return -1.0
}

// chooseBin picks which bin should be used for the FastRejectionSampler, where the choice of bin is weighted by its relative contribution to the overall integral of f.
func chooseBin(sumHeights float64, binHeights []float64, seed *rand.Rand) int {

	var cumulative float64 = 0.0
	for i := 0; i < len(binHeights); i++ {
		cumulative += binHeights[i] / sumHeights
		if cumulative > seed.Float64() {
			return i
		}
	}
	log.Fatalf("Error in chooseBin: failed to choose a bin.")
	return -1
}

// RejectionSample returns simulated values from an arbitrary function between a specified left and right bound using a simple rejection sampling method.
func RejectionSample(xLeft float64, xRight float64, yMax float64, f func(float64) float64, maxIteration int, seed *rand.Rand) float64 {
	var x, y float64
	for i := 0; i < maxIteration; i++ {
		x = RandFloat64InRange(xLeft, xRight, seed) //rand float64 in range xleft to xright
		y = f(x)
		if RandFloat64InRange(0.0, yMax, seed) < y {
			return y
		}
	}
	log.Fatalf("Exceeded max iterations.")
	return -1.0
}

// BoundedRejectionSample returns a rejection sample of a function f using a bounding function boundingSampler between a specified left and right bound.
func BoundedRejectionSample(boundingSampler func() (float64, float64), f func(float64) float64, xLeft float64, xRight float64, maxIteration int, seed *rand.Rand) (float64, float64) {
	var xSampler, ySampler, y float64
	for i := 0; i < maxIteration; i++ {
		xSampler, ySampler = boundingSampler()
		y = f(xSampler)
		if y > ySampler {
			log.Fatalf("BoundedRejectionSample: function was not a valid bounding function, ySampler is greater than y. xSampler: %e. ySampler: %e. y: %e.", xSampler, ySampler, y)
		}
		if RandFloat64InRange(0.0, ySampler, seed) < y {
			return xSampler, y
		}
	}
	log.Fatalf("BoundedRejectionSample: Exceeded max iteration.")
	return -1.0, -1.0
}

// ScaledBetaSampler returns an instatiation of RandBeta where the returned density has been scaled by the input variable 'multiplier'.
func ScaledBetaSampler(a float64, b float64, multiplier float64, seed *rand.Rand) func() (float64, float64) {
	return func() (float64, float64) {
		answer := RandBeta(a, b, seed)
		return answer, multiplier * BetaDist(answer, a, b)
	}
}

// BetaSampler returns an instantiation of RandBeta for a specified a and b parameter.
func BetaSampler(a float64, b float64, seed *rand.Rand) func() (float64, float64) {
	return func() (float64, float64) {
		answer := RandBeta(a, b, seed)
		return answer, BetaDist(answer, a, b)
	}
}

// RandGamma returns a random x,y point drawn from a gamma distribution with parameters alpha and beta. y corresponds to the function density at that x value.
// a > 1 uses the method from Marsaglia and Tsang 2000. Written for k, theta parameters, so the first step converts b to 1 / b to evaluate gamma in terms of alpha and beta parameters.
// a < 1 uses the method from Ahrens, J.H. and Dieter, U. (1974). Computer methods for sampling from gamma, beta, poisson and binomial distributions. Computing, 12, 223-246.
func RandGamma(a float64, b float64, seed *rand.Rand) (float64, float64) {
	if a < 0 || b < 0 {
		log.Fatalf("Error: The gamma distribution is defined with alpha and beta parameters greater than zero.")
	}
	b = 1 / b //convert parameter system
	var x, v, u, rExp float64
	if a < 1 {
		/* Marsaglia and Tsang code, does not appear to work.
		u = rand.Float64()
		return RandGamma(1.0+a, b) * math.Pow(u, 1.0/a)
		*/
		e1 := 0.36787944117144232159 //exp(-1), left as a constant to speed up computation
		e := 1.0 + e1*a
		for 1 > 0 { //repeat loop until breaks
			p := e * seed.Float64()
			rExp, _ = RandExp()
			if p >= 1.0 {
				x = -1 * math.Log((e-p)/a)
				if rExp >= (1.0-a)*math.Log(x) {
					break
				}
			} else {
				x = math.Exp(math.Log(p) / a)
				if rExp >= x {
					break
				}
			}
		}
		return b * x, GammaDist(a, b, b*x)
	}

	var d float64 = a - (1.0 / 3.0)
	var c float64 = (1.0 / 3.0) / math.Sqrt(d)
	for 1 > 0 {
		x = seed.NormFloat64()
		v = 1.0 + c*x
		for v <= 0 { //do while loop
			x = seed.NormFloat64()
			v = 1.0 + c*x
		}
		v = v * v * v
		u = seed.Float64()
		if u < 1-0.0331*x*x*x*x {
			break
		}
		if math.Log(u) < 0.5*x*x+d*(1-v+math.Log(v)) {
			break
		}
	}
	return b * d * v, GammaDist(a, b, b*d*v)
}

// GammaSampler returns an instantiation of RandGamma for specified a and b parameters.
func GammaSampler(a float64, b float64, seed *rand.Rand) func() (float64, float64) {
	return func() (float64, float64) {
		return RandGamma(a, b, seed)
	}
}
