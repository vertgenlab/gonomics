package numbers

import (
	"log"
	"math"
	"math/rand"
)

//returns a simulated value from a normal distribution.
func SampleInverseNormal(mu float64, sigma float64) float64 {
	return rand.NormFloat64()*sigma + mu
}

//pass in support of function
//truncate for infinite functions
func RejectionSample(xleft float64, xright float64, f func(float64) float64, maxIteration int) float64 {
	var x, y float64
	for i := 0; i < maxIteration; i++ {
		x = rand.Float64() //rand float64 in range xleft to xright
		y = f(x)
		if rand.Float64() < y {
			return y
		}
	}
	log.Fatalf("Exceeded max iterations.")
	return -1.0
}

/*
Returns a random number drawn from a gamma distribution with parameters alpha and beta.
Using the method from Marsaglia and Tsang 2000.
Written for k, theta parameters, so the first step converts b to 1 / b to evaluate gamma in terms of alpha and beta parameters.
*/
func RandGamma(a float64, b float64) float64 {
	if a < 0 || b < 0 {
		log.Fatalf("Error: The gamma distribution is defined with alpha and beta parameters greater than zero.")
	}
	b = 1 / b //convert parameter system
	var x, v, u float64
	if a < 1 {
		u = rand.Float64()
		return RandGamma(1.0+a, b) * math.Pow(u, 1.0/a)
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
	return b * d * v
}
