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
// TODO: memory allocation could be reduced since only the last row of R needs to be saved
func rombergsMethod(f func(float64) float64, a float64, b float64, precision float64, maxIter int) float64 {
	var n, m, i int
	var k float64

	var h []float64 = make([]float64, maxIter)
	var R [][]float64 = make([][]float64, maxIter)
	for i = 0; i < maxIter; i++ {
		R[i] = make([]float64, i+1)
	}

	R[0][0] = 0.5 * (f(a) + f(b))
	for n = 1; n < maxIter; n++ {
		h[n] = math.Exp2(float64(-n)) * (b - a)
		for k = 1; k <= math.Exp2(float64(n-1)); k++ {
			R[n][0] += f(a + (2*k-1)*h[n])
		}
		R[n][0] *= h[n]
		R[n][0] += 0.5 * R[n-1][0]
		for m = 1; m <= n; m++ {
			R[n][m] = R[n][m-1] + 1/(math.Pow(4, float64(m))-1)*(R[n][m-1]-R[n-1][m-1])
		}
		if R[n][n]-R[n][n-1] < precision {
			return R[n][n]
		}
	}
	log.Fatal("Error: Romberg's method did not converge.")
	return (0)
}

// DefiniteIntegral computes the definite integral of f(x) dx from start to end
func DefiniteIntegral(f func(float64) float64, start float64, end float64) float64 {
	return rombergsMethod(f, start, end, 1e-10, 20)
}
