package numbers

import (
//DEBUG: "fmt"
//DEBUG"log"
)

//BinomialDistLog returns log(BinomialDist), where log is the natural logarithm.
//This is ideal for very small probabilities to avoid underflow.
func BinomialDistLog(n int, k int, p float64) float64 {
	coefficient := BinomCoefficientLog(n, k)
	expression := BinomialExpressionLog(n, k, p)
	return MultiplyLog(coefficient, expression)
}

func BinomialExpressionLog(n int, k int, p float64) float64 {
	s := LogPow(p, float64(k))
	f := LogPow(1.0-p, float64(n-k))
	//DEBUG:log.Printf("n: %d. k: %d. p: %f. s: %e. f: %e.\n",n, k, p, s, f)
	return MultiplyLog(s, f)
}

//BinomialDistLogMap returns log(BinomialDist), where log is the natural logarithm.
//This function is similar to BinomialDistLog but passes in a map[int][]float64, where the int key
//refers to n and the []float64 map values are the corresponding binomial coefficients for index k.
//Useful to not recalculate the binomial coefficient each time when binomial densities must be constantly evaluated in logSpace, like in MCMC.
func BinomialDistLogSlice(n int, k int, p float64, binomCache [][]float64) float64 {
	//DEBUG: fmt.Printf("n: %v. k: %v. len(binomMap): %v.", n, k, len(binomMap))
	if binomCache[n] == nil {
		binomCache[n] = AddBinomMapEntry(n)
	}
	s := LogPow(p, float64(k))
	f := LogPow(1.0-p, float64(n-k))
	expression := MultiplyLog(s, f)
	return MultiplyLog(expression, binomCache[n][k])
}

//AddBinomMapEntry adds an entry to a binomMap containing a slice of binomial coefficients in logSpace for a particular n value.
func AddBinomMapEntry(n int) []float64 {
	var answer []float64
	answer = make([]float64, n+1)
	for k := 1; k < n+1; k++ {
		answer[k] = BinomCoefficientLog(n, k)
	}
	//DEBUG:fmt.Printf("Answer: %v\n", answer)
	return answer
}
