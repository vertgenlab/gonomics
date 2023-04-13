package numbers

import (
	"log"
	"math"

	"github.com/vertgenlab/gonomics/numbers/logspace"
)

// LogIntegrate evaluates log(int_a^b f(x)dx) in cases where f returns log(f(x)). Uses the rectangle rule.
func LogIntegrate(f func(float64) float64, a float64, b float64, n int) float64 {
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
	answer = logspace.Multiply(logspace.Average(f(currLeft), nextLeftEval), logDeltaX)
	var rightEval float64

	for i := 1; i < n; i++ {
		currLeft += deltaX
		currRight += deltaX
		rightEval = f(currRight)
		answer = logspace.Add(answer, logspace.Multiply(logspace.Average(nextLeftEval, rightEval), logDeltaX))
		nextLeftEval = rightEval
	}
	return answer
}

func LogIntegrateIterative(f func(float64) float64, a float64, b float64, maxIter int, relativeError float64) float64 {
	if maxIter < 2 {
		log.Fatalf("maxIterations for LogIntegrateIterative must be at least 2.")
	}
	if relativeError <= 0 {
		log.Fatalf("relativeError for LogIntegrateIterative must be greater than 0.")
	}
	n := 1000
	var prev, curr float64
	prev = LogIntegrate(f, a, b, n)

	for i := 0; i < maxIter; i++ {
		n := n * 10
		curr = LogIntegrate(f, a, b, n)
		if math.Abs(prev-curr)/curr < relativeError {
			//DEBUG: log.Printf("In LogIntegrateIterative: i=%v.", i)
			return curr
		}
		prev = curr
	}
	log.Fatalf("LogIntegrateIterative failed to converge below relative error: %f in maxIter: %v.", relativeError, maxIter)
	return (0)
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

// DefiniteSmallIntegral is like DefiniteIntegral with absolute error set to zero, so only relative error defines convergence conditions.
// slower than DefiniteIntegral, but more accurate for small values.
func DefiniteSmallIntegral(f func(float64) float64, start float64, end float64) float64 {
	return rombergsMethod(f, start, end, 0, 1e-6, 30)
}

// adaptiveSimponsHelper is the recursive core function for AdaptiveSimpsons
func adaptiveSimpsonsHelper(f func(float64) float64, a, b, midpoint, fa, fb, fMidpoint, wholeEstimate, errorThresh float64, maxDepth int) float64 {
	var h, leftMidpoint, rightMidpoint, fLeftMidpoint, fRightMidpoint, leftEstimate, rightEstimate, delta float64
	h = (b - a) / 2
	leftMidpoint = (a + midpoint) / 2
	rightMidpoint = (midpoint + b) / 2

	if maxDepth < 0 {
		log.Fatalf("Error in integration: exceeded maximum depth\n")
	} else if errorThresh/2 == errorThresh {
		log.Fatalf("Error in integration: the error threshold has gotten too small after many recursive calls\n")
	} else if a == leftMidpoint {
		log.Fatalf("Error in integration: the left side and midpoint have gotten too close to each other\n")
	}

	fLeftMidpoint = f(leftMidpoint)
	fRightMidpoint = f(rightMidpoint)
	leftEstimate = (h / 6) * (fa + 4*fLeftMidpoint + fMidpoint)
	rightEstimate = (h / 6) * (fMidpoint + 4*fRightMidpoint + fb)
	delta = leftEstimate + rightEstimate - wholeEstimate

	// Lyness 1969 + Richardson extrapolation; see article
	if math.Abs(delta) <= 15*errorThresh {
		return leftEstimate + rightEstimate + delta/15
	} else {
		return adaptiveSimpsonsHelper(f, a, midpoint, leftMidpoint, fa, fMidpoint, fLeftMidpoint, leftEstimate, errorThresh/2, maxDepth-1) + adaptiveSimpsonsHelper(f, midpoint, b, rightMidpoint, fMidpoint, fb, fRightMidpoint, rightEstimate, errorThresh/2, maxDepth-1)
	}
}

// AdaptiveSimpsons returns the integral from a to b of function f
// The error in the calculation should be less than or equal to errorThreshold.  If this can not be
// achieved within maxDepth number recursions, then the function aborts.
func AdaptiveSimpsons(f func(float64) float64, a float64, b float64, errorThreshold float64, maxDepth int) float64 {
	var midpoint, h, fa, fb, fMidpoint, s float64
	h = b - a
	midpoint = (a + b) / 2
	fa = f(a)
	fb = f(b)
	fMidpoint = f(midpoint)
	s = (h / 6) * (fa + 4*fMidpoint + fb)
	return adaptiveSimpsonsHelper(f, a, b, midpoint, fa, fb, fMidpoint, s, errorThreshold, maxDepth)
}

// adaptiveSimponsLogHelper is the recursive core function for AdaptiveSimpsonsLog
func adaptiveSimpsonsLogHelper(f func(float64) float64, a, b, midpoint, fa, fb, fMidpoint, wholeEstimate, errorThresh float64, maxDepth int) float64 {
	const logFour float64 = 1.386294
	const logFifteen float64 = 2.70805
	const logHalf float64 = -0.6931472
	var estimateFromHalves, logHOverSix, h, leftMidpoint, rightMidpoint, fLeftMidpoint, fRightMidpoint, leftEstimate, rightEstimate, delta float64
	h = (b - a) / 2
	leftMidpoint = (a + midpoint) / 2
	rightMidpoint = (midpoint + b) / 2

	if maxDepth < 0 {
		log.Fatalf("Error in integration: exceeded maximum depth\n")
	} else if logspace.Multiply(errorThresh, logHalf) == errorThresh {
		log.Fatalf("Error in integration: the error threshold has gotten too small after many recursive calls\n")
	} else if a == leftMidpoint {
		log.Fatalf("Error in integration: the left side and midpoint have gotten too close to each other. a: %e. b: %e. Midpoint: %e. LeftMidpoint:%e. Fa: %e. Fb: %e. MaxDepth: %d.\n", a, b, midpoint, leftMidpoint, fa, fb, maxDepth)
	}

	fLeftMidpoint = f(leftMidpoint)
	fRightMidpoint = f(rightMidpoint)
	//DEBUG: log.Printf("fLeftMidpoint: %e. fRightMidpoint: %e.\n", fLeftMidpoint, fRightMidpoint)
	logHOverSix = math.Log(h / 6)
	leftEstimate = logspace.Multiply(logHOverSix, logspace.Add(logspace.Add(fa, logspace.Multiply(logFour, fLeftMidpoint)), fMidpoint))
	rightEstimate = logspace.Multiply(logHOverSix, logspace.Add(logspace.Add(fMidpoint, logspace.Multiply(logFour, fRightMidpoint)), fb))
	estimateFromHalves = logspace.Add(leftEstimate, rightEstimate)

	//log.Printf("maxDepth:%d, left:%f, right:%f, fromHalves:%f, whole:%f\n", maxDepth, math.Exp(leftEstimate), math.Exp(rightEstimate), math.Exp(estimateFromHalves), math.Exp(estimateFromHalves))

	if estimateFromHalves > wholeEstimate {
		delta = logspace.Subtract(estimateFromHalves, wholeEstimate)
		if delta <= logspace.Multiply(logFifteen, errorThresh) {
			return logspace.Add(logspace.Add(leftEstimate, rightEstimate), logspace.Divide(delta, logFifteen))
		}
	} else if wholeEstimate > estimateFromHalves {
		delta = logspace.Subtract(wholeEstimate, estimateFromHalves)
		if delta <= logspace.Multiply(logFifteen, errorThresh) {
			return logspace.Add(logspace.Add(leftEstimate, rightEstimate), logspace.Divide(delta, logFifteen))
		}
	}
	return logspace.Add(adaptiveSimpsonsLogHelper(f, a, midpoint, leftMidpoint, fa, fMidpoint, fLeftMidpoint, leftEstimate, logspace.Multiply(errorThresh, logHalf), maxDepth-1), adaptiveSimpsonsLogHelper(f, midpoint, b, rightMidpoint, fMidpoint, fb, fRightMidpoint, rightEstimate, logspace.Multiply(errorThresh, logHalf), maxDepth-1))
}

// AdaptiveSimpsons returns the log of the integral from a to b of g(x), where f(x) = log(g(x))
// The error in the calculation should be less than or equal to errorThreshold.  If this can not be
// achieved within maxDepth number recursions, then the function aborts.
func AdaptiveSimpsonsLog(f func(float64) float64, a float64, b float64, errorThreshold float64, maxDepth int) float64 {
	const logFour float64 = 1.386294
	var midpoint, h, fa, fb, fMidpoint, s float64
	h = b - a
	midpoint = (a + b) / 2
	fa = f(a)
	fb = f(b)
	fMidpoint = f(midpoint)
	s = logspace.Multiply(math.Log(h/6), logspace.Add(logspace.Add(fa, logspace.Multiply(logFour, fMidpoint)), fb))
	return adaptiveSimpsonsLogHelper(f, a, b, midpoint, fa, fb, fMidpoint, s, math.Log(errorThreshold), maxDepth)
}
