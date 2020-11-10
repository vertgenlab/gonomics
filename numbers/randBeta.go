package numbers

import (
	"math"
	"math/rand"
	"log"
)

const RelError float64 = 1.0e-8
const Small float64 = 1.0e-30
const MaxIterations int = 200

//RandBeta returns a random number drawn from a beta distribution with parameters alpha and beta as an xy point where y is the function density.
//Implemented from principles in R. C. H. Cheng (1978). Generating beta variates with nonintegral shape parameters. Communications of the ACM 21, 317-322.
func RandBeta(alpha float64, beta float64) float64 {
	//special returns
	if alpha < 0 || beta < 0 {
		log.Fatalf("Error in RandBeta: alpha and beta must be greater than 0.  alpha: %f. beta: %f.", alpha, beta)
	}
	if math.IsInf(alpha, 1) && math.IsInf(beta, 1) {
		return 0.5
	}
	if alpha == 0 && beta == 0 {
		return 1.0
	}
	if math.IsInf(alpha, 1) || beta == 0 {
		return 1.0
	}
	if math.IsInf(beta, 1) || alpha == 0 {
		return 0.0
	}

	var sum, a, b, bet, delta, k1, k2, r, s, t, u1, u2, v, w, y, z float64
	var i int
	a = MinFloat64(alpha, beta)
	b = MaxFloat64(alpha, beta)
	sum = a + b

	if a <= 1.0 {
		//initializations
		bet = 1.0 / a
		delta = 1.0 + b - a
		k1 = delta * (0.0138889 + 0.0416667 * a) / (b * bet - 0.777778)
		k2 = 0.25 + (0.5 + 0.25 / delta) * a

		for i = 0; i < MaxIterations; i++ {
			u1 = rand.Float64()
			u2 = rand.Float64()
			if u1 < 0.5 {
				y = u1 * u2
				z = u1 * y
				if (0.25 * u2 + z - y >= k1) {
					continue
				}
			} else {
				z = u1 * u1 * u2
				if (z <= 0.25) {
					//this block of code was broken into a helper function in the R implementation
					v = beta * math.Log(u1 / (1.0 - u1))
					if (v <= math.Log(math.MaxFloat64)) {
						w = alpha * math.Exp(v)
						if math.IsInf(w, 1) {
							w = math.MaxFloat64
						}
					} else {
						w = math.MaxFloat64
					}
					break
				}
				if z >= k2 {
					continue
				}
			}

			v = beta * math.Log(u1 / (1.0 - u1))
			if (v <= math.Log(math.MaxFloat64)) {
				w = b * math.Exp(v)
				if math.IsInf(w, 1) {
					w = math.MaxFloat64
				}
			} else {
				w = math.MaxFloat64
			}

			if (sum * (math.Log(sum / (a + w)) + v) - 1.3862944 >= math.Log(z)) {
				break
			}
		}

		if i == MaxIterations {
			log.Fatalf("Failed to converge.")
			return -1.0
		}

		if alpha == a {
			return a / (a + w)
		}
		return w / (a + w)
	} else {//else we use algorithm BB from page 320.
		var u1, u2 float64
		var bet float64 = math.Sqrt((sum - 2.0) / (2.0 * a * b - sum))
		var gamma float64 = a + 1.0 / bet
		var firstTime bool = true
		for firstTime || r + sum * math.Log(sum / (b + w)) < t {//dowhile
			firstTime = false
			u1 = rand.Float64()
			u2 = rand.Float64()

			v = beta * math.Log(u1 / (1.0 - u1))
			if (v <= math.Log(math.MaxFloat64)) {
				w = a * math.Exp(v)
				if math.IsInf(w, 1) {
					w = math.MaxFloat64
				}
			} else {
				w = math.MaxFloat64
			}

			z = u1 * u1 * u2
			r = gamma * v - 1.3862944
			s = a + r - w
			if (s + 2.609438 >= 5.0 * z) {
				break
			}
			t = math.Log(z)
			if s > t {
				break
			}
		}
		if alpha != a {
			return b / (b + w)
		} else {
			return w / (b + w)
		}
	}
}

//Translation of an implementation of the Regularized Incomplete Beta function Copyright (c) 2016, 2017 Lewis Van Winkle, zlib license. https://codeplea.com/incomplete-beta-function-c
func incompleteBetaHelper(a float64, b float64, x float64) float64 {
	if (x > (a+1.0)/(a+b+2.0)) {
		return (1.0-incompleteBetaHelper(b, a, 1.0-x))
	}
	logBeta := math.Log(BetaFunc(a, b))
	front := math.Exp(math.Log(x)*a + math.Log(1.0-x)*b-logBeta) / a

	var f, c, d float64 = 1.0, 1.0, 0.0

	var numerator, m float64
	var i int
	for i = 0; i <= MaxIterations; i++ {
		m = float64(i / 2)
		if i == 0 {
			numerator = 1.0
		} else if (i % 2 == 0) {
			numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m))
		} else {
			numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1))
		}

		d = 1.0 + numerator * d
		if (math.Abs(d) < Small) {
			d = Small
		}

		d = 1.0 / d
		c = 1.0 + numerator / c

		if (math.Abs(c) < Small) {
			c = Small
		}

		f *= c * d

		if math.Abs(1.0-(c*d)) < RelError {
			return front * (f-1.0)
		}
	}
	log.Fatalf("Failed to converge.")
	return -1.0
}