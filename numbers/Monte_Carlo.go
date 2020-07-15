package numbers

import (
	"math/rand"
	//"math"
	"log"
)

//returns a simulated value from a normal distribution.
func SampleInverseNormal(mu float64, sigma float64) float64 {
	return rand.NormFloat64()*sigma + mu
}

//pass in support of function
//truncate for infinite functions
func RejectionSample(xleft float64, xright float64, f func(float64) float64, maxIteration int) float64 {
	var x, y, yRand float64
	for i := 0; i < maxIteration; i++ {
		x = rand.Float64() //rand float64 in range xleft to xright
		y = f(x)
		yRand = rand.Float64() //accept reject random variable
		if yRand < y {
			return y
		}
	}
	log.Fatalf("Exceeded max iterations.")
	return -1.0
}
