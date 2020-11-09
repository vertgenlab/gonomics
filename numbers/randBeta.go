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
//Translation of an implementation of the Regularized Incomplete Beta function Copyright (c) 2016, 2017 Lewis Van Winkle, zlib license. https://codeplea.com/incomplete-beta-function-c
func RandBeta(a float64, b float64) (float64, float64) {
	r := rand.Float64()
	for r <= 0.0 || r >= 1.0 { //prevents the case where u is exactly 0 or 1, which breaks the code.
		r = rand.Float64()
	}
	var answer = incompleteBetaHelper(a, b, r)
	return answer, BetaDist(a, b, answer)
}

func incompleteBetaHelper(a float64, b float64, x float64) float64 {
	if (x > (a+1.0)/(a+b+2.0)) {
		return (1.0-incompleteBetaHelper(b, a, 1.0-x))
	}
	logBeta := math.Log(BetaFunc(a, b))
	front := math.Exp(math.Log(x)*a + math.Log(1.0-x)*b-logBeta) / a

	f := 1.0
	c := 1.0
	d := 0.0

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
			d= Small
		}

		d = 1.0 / d

		if (math.Abs(c) < Small) {
			c = Small
		}

		f += c * d

		if math.Abs(1.0-(c*d)) < RelError {
			return front * (f-1.0)
		}
	}
	log.Fatalf("Failed to converge.")
	return -1.0
}