// Package numbers provides a variety of functions related to statistics
package numbers

import (
	"log"
	"math"

	"github.com/vertgenlab/gonomics/numbers/logspace"
)

// carefulMultDivFloat tries to gracefully handle the case when you have a
// bunch of numbers being multiplied together in a numerator (numer)
// and then being divided by a bunch of numbers
// in a denominator (denom).  For example numer={2,3,4} and denom={6,3,4}
// would be equal to 2*3*4/6/3/4.
func carefulMultDivFloat(numer []float64, denom []float64) float64 {
	var answer float64 = 1
	var i, j int = 0, 0

	for i < len(numer) || j < len(denom) {
		if (answer <= 1 && i < len(numer)) || j == len(denom) {
			if math.MaxFloat64/numer[i] < answer {
				log.Fatal("Error: carefulMultDiv detected that it would have had overflow\n")
			}
			answer = answer * numer[i]
			i++
		} else {
			if math.SmallestNonzeroFloat64*denom[j] > answer {
				if i == len(numer) {
					return math.SmallestNonzeroFloat64
				} else {
					log.Fatal("Error: carefulMultDiv detected that it would have had underflow\n")
				}
			}
			answer = answer / denom[j]
			j++
		}
	}
	return answer
}

// carefulMultDivInt does the same thing as carefulMultDivFloat, but for
// ints and raises a fatal error if the answer can not be represented as an int
// TODO: right now the overflow check assumes we are on 64bit arch.
func carefulMultDivInt(numer []int, denom []int) int {
	var answer int = 1
	var i, j int = 0, 0

	for i < len(numer) || j < len(denom) {
		// if there are still numbers in the denominator
		// and the number in the denominator goes in evenly to
		// the running calculation
		if j < len(denom) && answer%denom[j] == 0 {
			answer = answer / denom[j]
			j++
		} else if i < len(numer) {
			if math.MaxInt64/numer[i] < answer {
				log.Fatal("Error: integer overflow was detected\n")
			} else {
				answer = answer * numer[i]
				i++
			}
		} else {
			log.Fatal("Error: integer division could not be done successfully\n")
		}
	}
	return answer
}

// matrix is in the form of:
// [a, b]
// [c, d]
// pvalue returned is for a being smaller than expected
// in relation to b, given the
// ratio of c to d.
func fisherExactLess(a, b, c, d int) float64 {
	var currProb float64 = fisherProbLess(a, b, c, d)
	var runningTotal float64 = currProb
	for a > 0 && d > 0 {
		a = a - 1
		b = b + 1
		c = c + 1
		d = d - 1
		//currProb = fisherProbLess(a, b, c, d) // this way may be more resistant to propagation of error, but slower
		currProb = currProb * float64(a+1) / float64(c) * float64(d+1) / float64(b)
		runningTotal += currProb
	}
	return runningTotal
}

func fisherProbLess(a, b, c, d int) float64 {
	var n int = a + b + c + d
	numer := make([]float64, n)
	denom := make([]float64, n)
	var i int = 0
	for w := a + 1; w < a+b+1; w++ {
		numer[i] = float64(w)
		i++
	}
	for x := d + 1; x < c+d+1; x++ {
		numer[i] = float64(x)
		i++
	}
	for y := c + 1; y < a+c+1; y++ {
		numer[i] = float64(y)
		i++
	}
	for z := b + 1; z < b+d+1; z++ {
		numer[i] = float64(z)
		i++
	}
	for j := 1; j < n+1; j++ {
		denom[j-1] = float64(j)
	}
	return carefulMultDivFloat(numer, denom)
}

// FisherExact computes a one-sided Fisher's
// Exact test on the 2x2 table provided
// The test is for the matrix:
// [a b]
// [c d]
// aSmall being true tests for the ratio of a to b
// being small, given the ratio of c to d
// aSmall being false tests for the ratio of a to b
// being large, given the ratio of c to d.
func FisherExact(a, b, c, d int, aSmall bool) float64 {
	if aSmall {
		return fisherExactLess(a, b, c, d)
	} else {
		return fisherExactLess(c, d, a, b)
	}
}

// BinomCoefficient calculates the Binomial Coefficient, which
// is also called the "Choose" Function. The answer
// returned is "n choose k" or n!/(n-k)!k!
func BinomCoefficient(n int, k int) int {
	if n < 0 || k < 0 || k > n {
		log.Fatalf("The binomial coefficient call could not be handled: n=%d and k=%d\n", n, k)
	}
	if n-k > k {
		k = n - k
	}
	// this special case is handled here so that we don't ask for negative memory for denom
	if k == n {
		return 1
	}
	numer := make([]int, 0, n-k)
	denom := make([]int, 0, n-k-1)
	var x, y int
	for x = k + 1; x < n+1; x++ {
		numer = append(numer, x)
	}
	for y = 2; y < n-k+1; y++ {
		denom = append(denom, y)
	}
	return carefulMultDivInt(numer, denom)
}

// BinomCoefficientLog returns log(n choose k), where log is the natural logarithm.
// Ideal for large numbers as this raises the overflow ceiling considerably.
func BinomCoefficientLog(n int, k int) float64 {
	if n < 0 || k < 0 || k > n {
		log.Fatalf("The binomial coefficient call could not be handled: n=%d and k=%d\n", n, k)
	}
	if n-k > k {
		k = n - k
	}
	// this special case is handled here so that we don't ask for negative memory for denom
	if k == n {
		return 0.0
	}
	var x, y int
	var numer, denom float64 = 0.0, 0.0
	for x = k + 1; x < n+1; x++ {
		numer = logspace.Multiply(numer, math.Log(float64(x)))
	}
	for y = 2; y < n-k+1; y++ {
		denom = logspace.Multiply(denom, math.Log(float64(y)))
	}
	return logspace.Divide(numer, denom)
}

// Factorial returns n! in normal space.
func Factorial(n int) int {
	return int(math.Gamma(float64(n + 1)))
}

// DigitsBaseTen returns the.
func DigitsBaseTen(x int) int {
	var count int = 1
	if x < 0 {
		x = -1 * x
		count++
	}
	for x >= 10 {
		x = x / 10
		count++
	}
	return count
}

// AbsInt returns the absolute value of an input of type int.
func AbsInt(x int) int {
	if x < 0 {
		return -x
	} else {
		return x
	}
}
