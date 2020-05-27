package numbers

import (
	"log"
	"math"
)

// There are a number of ways to evaluate a definite integral computationally.
// Romberg's method seems like a good mix of accuracy and coding difficulty,
// but there are better methods out there if more speed or accuracy are needed.
// This code tries to follow the algorithm and variable names used here:
// https://en.wikipedia.org/wiki/Romberg's_method
func rombergsMethod(f func(float64) float64, a float64, b float64, estimatedError float64, maxIter int) float64 {
	var n, m int
	var kMax, k float64

	var h float64
	var currR, prevR []float64 = make([]float64, maxIter), make([]float64, maxIter)

	prevR[0] = 0.5 * (f(a) + f(b))
	for n = 1; n < maxIter; n++ {
		// compute the current h value
		h = math.Exp2(float64(-n)) * (b - a)

		// compute R[n][0]
		currR[0] = 0 // needed because memory is being reused
		kMax = math.Exp2(float64(n - 1))
		for k = 1; k <= kMax; k++ {
			currR[0] += f(a + (2*k-1)*h)
		}
		currR[0] *= h
		currR[0] += 0.5 * prevR[0]

		// now that we have R[n][0], we can compute R[n][m] where m > 0
		for m = 1; m <= n; m++ {
			currR[m] = currR[m-1] + 1/(math.Pow(4, float64(m))-1)*(currR[m-1]-prevR[m-1])
		}

		// now checking to see if we have convergence
		// some people use R[n][n]-R[n][n-1]
		// and some use R[n][n]-R[n-1][n-1]
		// these appear to be related by a constant of 1/(4^n-1) with
		// R[n][n]-R[n-1][n-1] being more conservative, so we will use that one
		if math.Abs(currR[n]-prevR[n-1]) < estimatedError {
			return currR[n]
		}

		// swap prev and curr so that current becomes prev
		prevR, currR = currR, prevR
	}
	log.Fatal("Error: Romberg's method did not converge.")
	return (0)
}

// DefiniteIntegral computes the definite integral of f(x) dx from start to end
func DefiniteIntegral(f func(float64) float64, start float64, end float64) float64 {
	return rombergsMethod(f, start, end, 1e-5, 30)
}
