package numbers

import (
	"log"
	"math"
	"math/rand"
	//DEBUG: "fmt"
)

//Sample InverseNormal returns a simulated value from a normal distribution.
func SampleInverseNormal(mu float64, sigma float64) float64 {
	return rand.NormFloat64()*sigma + mu
}

//InitializeFastRejectionSampler takes in the parameters of a rejection sampler and returns the binHeights and sumHeights variables.
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
			binHeights[i] = MaxFloat64(fCurrLeft, fCurrRight)
		} else {
			fCurrLeft = fCurrRight
			currRight += stepSize
			fCurrRight = f(currRight)
			binHeights[i] = MaxFloat64(fCurrLeft, fCurrRight)
		}
		sumHeights += binHeights[i]
	}
	return binHeights, sumHeights
}

//FastRejectionSampler returns simulated values from an a func(float64) float64 between a left and right value using an optimized rejection sampler
//that divides the function support into discrete bins with optimized sampling heights.
//maxSampleDepth triggers the log.Fatalf in the RejectionSample func, and samples is the number of values to be returned.
func FastRejectionSampler(xLeft float64, xRight float64, f func(float64) float64, bins int, maxSampleDepth int, samples int) []float64 {
	var answer []float64 = make([]float64, samples)
	var stepSize float64 = (xRight - xLeft) / float64(bins)
	binHeights, sumHeights := InitializeFastRejectionSampler(xLeft, xRight, f, bins)
	for j := 0; j < samples; j++ {
		answer[j] = RejectionSampleChooseBin(xLeft, xRight, stepSize, f, maxSampleDepth, sumHeights, binHeights)
	}
	return answer
}

func RejectionSampleChooseBin(xLeft float64, xRight float64, stepSize float64, f func(float64) float64, maxIteration int, sumHeights float64, binHeights []float64) float64 {
	var x, y float64
	var currBin int
	var currLeft, currRight float64
	for i := 0; i < maxIteration; i++ {
		currBin = chooseBin(sumHeights, binHeights)
		currLeft = xLeft + float64(currBin)*stepSize
		currRight = currLeft + stepSize
		x = RandFloat64InRange(currLeft, currRight)
		y = f(x)
		if RandFloat64InRange(0.0, binHeights[currBin]) < y {
			return x
		}
	}
	log.Fatalf("Exceeded max iteration in RejectionSampleChooseBin.")
	return -1.0
}

//chooseBin picks which bin should be used for the FastRejectionSampler, where the choice of bin is weighted by its relative contribution to the overall integral of f.
func chooseBin(sumHeights float64, binHeights []float64) int {
	var rand float64 = rand.Float64()
	var cumulative float64 = 0.0
	for i := 0; i < len(binHeights); i++ {
		cumulative += (binHeights[i] / sumHeights)
		if cumulative > rand {
			return i
		}
	}
	log.Fatalf("Error in chooseBin: failed to choose a bin.")
	return -1
}

//RejectionSample returns simulated values from an arbitrary function between a specified left and right bound using a simple rejection sampling method.
func RejectionSample(xLeft float64, xRight float64, yMax float64, f func(float64) float64, maxIteration int) float64 {
	var x, y float64
	for i := 0; i < maxIteration; i++ {
		x = RandFloat64InRange(xLeft, xRight) //rand float64 in range xleft to xright
		y = f(x)
		if RandFloat64InRange(0.0, yMax) < y {
			return y
		}
	}
	log.Fatalf("Exceeded max iterations.")
	return -1.0
}

func BoundedRejectionSample(boundingSampler func() (float64, float64), f func(float64) float64, xLeft float64, xRight float64, maxIteration int) float64 {
	var xSampler, ySampler, y float64
	for i := 0; i < maxIteration; i++ {
		xSampler, ySampler = boundingSampler()
		y = f(xSampler)
		if y > ySampler {
			log.Fatalf("BoundedRejectionSample: function was not a valid bounding function, ySampler is greater than y.")
		}
		if RandFloat64InRange(0.0, ySampler) < y {
			return y
		}
	}
	log.Fatalf("BoundedRejectionSample: Exceeded max iteration.")
	return -1.0
}

//RandExp Returns a random variable as a float64 from a standard exponential distribution. f(x)=e**-x.
//Algorithm from Ahrens, J.H. and Dieter, U. (1972). Computer methods for sampling from the exponential and normal distributions. Comm. ACM, 15, 873-882.
func RandExp() (float64, float64) {
	//q series where q[k-1] = sum(log(2)^k / k!) for k=1,2,...n
	q := [16]float64{0.6931471805599453, 0.9333736875190459, 0.9888777961838675, 0.9984959252914960, 0.9998292811061389, 0.9999833164100727, 0.9999985691438767, 0.9999998906925558, 0.9999999924734159, 0.9999999995283275, 0.9999999999728814, 0.9999999999985598, 0.9999999999999289, 0.9999999999999968, 0.9999999999999999, 1.0000000000000000}

	var a float64 = 0.0
	var r float64 = rand.Float64()
	for r <= 0.0 || r >= 1.0 { //prevents the case where u is exactly 0 or 1, which breaks the code.
		r = rand.Float64()
	}

	for 1 > 0 {
		r += r
		if r > 1.0 {
			break
		}
		a += q[0]
	}
	r -= 1
	if r <= q[0] {
		return a + r, ExpDist(a + r)
	}

	var i int = 0
	ustart := rand.Float64()
	umin := ustart

	for r > q[i] {
		ustart = rand.Float64()
		if umin > ustart {
			umin = ustart
		}
		i++
	}
	return a + umin*q[0], ExpDist(a + umin*q[0])
}

func BetaSampler(a float64, b float64) func() (float64, float64) {
	return func() (float64, float64) {
		answer := RandBeta(a, b)
		return answer, BetaDist(answer, a, b)
	}
}

//RandGamma returns a random x,y point drawn from a gamma distribution with parameters alpha and beta. y corresponds to the function density at that x value.
//a > 1 uses the method from Marsaglia and Tsang 2000. Written for k, theta parameters, so the first step converts b to 1 / b to evaluate gamma in terms of alpha and beta parameters.
//a < 1 uses the method from Ahrens, J.H. and Dieter, U. (1974). Computer methods for sampling from gamma, beta, poisson and binomial distributions. Computing, 12, 223-246.
func RandGamma(a float64, b float64) (float64, float64) {
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
			p := e * rand.Float64()
			rExp, _ = RandExp()
			if p >= 1.0 {
				x = -1 * math.Log((e-p)/a)
				if rExp >= (1.0-a)*math.Log(x) {
					break
				}
			} else {
				x = math.Exp(math.Log(p) / a)
				rExp, _ = RandExp()
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
		x = rand.NormFloat64()
		v = 1.0 + c*x
		for v <= 0 { //do while loop
			x = rand.NormFloat64()
			v = 1.0 + c*x
		}
		v = v * v * v
		u = rand.Float64()
		if u < 1-0.0331*x*x*x*x {
			break
		}
		if math.Log(u) < 0.5*x*x+d*(1-v+math.Log(v)) {
			break
		}
	}
	return b * d * v, GammaDist(a, b, b*d*v)
}

func GammaSampler(a float64, b float64) func() (float64, float64) {
	return func() (float64, float64) {
		return RandGamma(a, b)
	}
}
