package numbers

/*
	This file contains an implementation of an algorithm described in pseudocode in the following academic publication.
	Generating Beta Variates with Nonintegral Shape Parameters.
	R. C. H. Cheng
	Communications of the ACM, April 1978.
*/

import (
	"math"
	"math/rand"
	"log"
)

const RelError float64 = 1.0e-8
const Small float64 = 1.0e-30
const MaxIterations int = 200
const Ln4 float64 = 1.38629436112

//RandBeta is the basic beta variate generator from Cheng 1978. Uses the BA algorithm, which is less optimized, but still runs effectively for my use case. 
//More optimized algorithms (with greater programming complexity) are described in the Cheng paper, and should be implemented (TODO) if required.
func RandBeta(a float64, b float64) float64 {
	var alpha float64 = a + b
	var beta, gamma, u1, u2, w, v float64
	if MinFloat64(a, b) <= 1 {
		beta = MaxFloat64(1.0 / a, 1.0 / b)
	} else {
		beta = math.Sqrt((alpha - 2.0)/(2 * a * b - alpha))
	}
	gamma = a + 1.0 / beta

	for i := 0; i < MaxIterations; i++ {
		u1 = rand.Float64()
		u2 = rand.Float64()
		v = beta * math.Log(u1 / (1 - u1))
		w = a * math.Exp(v)

		if alpha * math.Log(alpha / (b + w)) + gamma * v - Ln4 < math.Log(u1 * u1 * u2) {
			continue
		}
		return w / (b + w)
	}
	log.Fatalf("RandBetaBA: Failed to find an accepted value in the max number of iterations.")
	return -1.0
}
