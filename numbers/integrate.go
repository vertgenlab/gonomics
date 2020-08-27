package numbers

import (
	"log"
	"math"
	//"fmt"
)

//logIntegrate evaluates log(int_a^b f(x)dx) in cases where f returns log(f(x)). Uses the rectangle rule.
func logIntegrate(f func(float64) float64, a float64, b float64, n int) float64 {
	if a >= b {
		log.Fatalf("logIntegrate failed, left bound must be smaller than right bound.")
	}
	var deltaX float64 = (b - a) / float64(n)
	var logDeltaX float64 = math.Log(deltaX)
	var currLeft float64 = a //this variable stores the left bound of the current rectangle.
	var currRight float64 = a + deltaX
	var answer float64
	//first time, sets answer as the area of the first rectangle
	var nextLeftEval float64 = f(currRight)
	answer = MultiplyLog(MidpointLog(f(currLeft), nextLeftEval), logDeltaX)
	var rightEval float64
	
	for i := 1; i < n; i++ {
		currLeft += deltaX
		currRight += deltaX
		rightEval = f(currRight)
		answer += MultiplyLog(MidpointLog(nextLeftEval, rightEval), logDeltaX)
		nextLeftEval = rightEval
	}
	return answer
}

// There are a number of ways to evaluate a definite integral computationally.
// Romberg's method seems like a good mix of accuracy and coding difficulty,
// but there are better methods out there if more speed or accuracy are needed.
// This code tries to follow the algorithm and variable names used here:
// https://en.wikipedia.org/wiki/Romberg's_method
func rombergsMethod(f func(float64) float64, a float64, b float64, estimatedError float64, relativeEstError float64, maxIter int) float64 {
	var n, m int
	var kMax, k, h, currEstError float64
	var currR, prevR []float64 = make([]float64, maxIter), make([]float64, maxIter)
	var minIter int = 10

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
		// log.Printf("prevEst=%e, currEst=%e\n", prevR[n-1], currR[n])
		currEstError = math.Abs(currR[n] - prevR[n-1])
		//fmt.Printf("currValue: %e. currError: %e\n", currR[n], currEstError)
		if (currEstError < estimatedError || currEstError < relativeEstError*math.Abs(currR[n])) && n >= minIter {
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
	return rombergsMethod(f, start, end, 1e-8, 1e-8, 30)
}

//DefiniteIntegral with absolute error set to zero, so only relative error defines convergence conditions.
//slower than DefiniteIntegral, but more accurate for small values.
func DefiniteSmallIntegral(f func(float64) float64, start float64, end float64) float64 {
	return rombergsMethod(f, start, end, 0, 1e-6, 30)
}
