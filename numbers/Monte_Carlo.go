package numbers

/*
	Implementation ported from the R source code. 
	https://github.com/wch/r-source/blob/trunk/src/nmath/rgamma.c
	https://github.com/wch/r-source/blob/trunk/src/nmath/sexp.c
*/

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

/* Returns a random variable as a float64 from a standard exponential distribution. f(x)=e**-x.
Algorithm from Ahrens, J.H. and Dieter, U. (1972). Computer methods for sampling from the exponential and normal distributions. Comm. ACM, 15, 873-882.
*/

func RandExp() float64 {
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
	if (r <= q[0]) {
		return a + r
	}

	var i int = 0
	ustart := rand.Float64()
	umin := ustart

	for (r > q[i]) {
		ustart = rand.Float64()
		if umin > ustart {
			umin = ustart
		}
		i++
	}
	return a + umin * q[0]
}

/*
Returns a random number drawn from a gamma distribution with parameters alpha and beta.
a > 1 uses the method from Marsaglia and Tsang 2000. Written for k, theta parameters, so the first step converts b to 1 / b to evaluate gamma in terms of alpha and beta parameters.
a < 1 uses the method from Ahrens, J.H. and Dieter, U. (1974). Computer methods for sampling from gamma, beta, poisson and binomial distributions. Computing, 12, 223-246.
*/

func RandGamma(a float64, b float64) float64 {
	if a < 0 || b < 0 {
		log.Fatalf("Error: The gamma distribution is defined with alpha and beta parameters greater than zero.")
	}
	b = 1 / b //convert parameter system
	var x, v, u float64
	if a < 1 {
		/* Marsaglia and Tsang code, does not appear to work.
		u = rand.Float64()
		return RandGamma(1.0+a, b) * math.Pow(u, 1.0/a)
		*/
		e1 := 0.36787944117144232159 //exp(-1), left as a constant to speed up computation
		e := 1.0 + e1 * a
		for 1 > 0 { //repeat loop until breaks
			p := e * rand.Float64()
			if p >= 1.0 {
				x = -1 * math.Log((e - p) / a)
				if RandExp() >= (1.0 - a) * math.Log(x) {
					break
				} 
			} else {
				x = math.Exp(math.Log(p) / a)
				if RandExp() >= x {
					break
				}
			}
		}
		return b * x
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
