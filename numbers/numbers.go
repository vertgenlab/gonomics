package numbers

import (
	"log"
	"math"
)

func carefulMultDiv(numer []float64, denom []float64) float64 {
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

// matrix is in the form of:
// [a, b]
// [c, d]
// pvalue returned is for a being smaller than expected
// in relation to b, given the
// ratio of c to d
func fisherExactLess(a, b, c, d int) float64 {
	var currProb float64 = fisherProbLess(a, b, c, d)
	var runningTotal float64 = currProb
	for a > 0 && d > 0 {
		a = a - 1
		b = b + 1
		c = c + 1
		d = d - 1
		//currProb = fisherProbLess(a, b, c, d) // this way may be more resistant to propogation of error, but slower
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
	return carefulMultDiv(numer, denom)
}

// test is for the matrix:
// [a b]
// [c d]
// aSmall being true tests for the ratio of a to b
// being small, given the ratio of c to d
// aSmall being false tests for the ratio of a to b
// being large, given the ratio of c to d
func FisherExact(a, b, c, d int, aSmall bool) float64 {
	if aSmall {
		return fisherExactLess(a, b, c, d)
	} else {
		return fisherExactLess(c, d, a, b)
	}
}
