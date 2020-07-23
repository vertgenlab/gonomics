package numbers

import (
	"github.com/vertgenlab/gonomics/common"
	"log"
	"math"
)

//TODO: Test functions for distribution values.

//returns the normal distribution value x for a distribution with mean mu and standard deviation sigma.
func NormalDist(x float64, mu float64, sigma float64) float64 {
	return (1 / (sigma * math.Sqrt(2*math.Pi))) * math.Exp(-0.5*math.Pow((x-mu)/sigma, 2))
}

func StandardNormalDist(x float64) float64 {
	return NormalDist(x, 0, 1)
}

func BinomialDist(n int, k int, p float64) float64 {
	return float64(BinomCoefficient(n, k)) * math.Pow(p, float64(k)) * math.Pow(1-p, float64(n-k))
}

func PoissonDist(k int, lambda float64) float64 {
	if k < 0 {
		log.Fatalf("The poisson distribution is supported for k > 0.")
	}
	if lambda <= 0 {
		log.Fatalf("The poisson distribution is supported for lambda >= 0.")
	}
	return (math.Pow(lambda, float64(k)) * math.Pow(math.E, -lambda)) / float64(Factorial(k))
}

func BetaDist(x float64, alpha float64, beta float64) float64 {
	if alpha <= 0 {
		log.Fatalf("Alpha parameter must be greater than 0.")
	}
	if beta <= 0 {
		log.Fatalf("Beta parameter must be greater than 0.")
	}
	if x < 0 || x > 1 {
		log.Fatalf("Value x out of range. The beta distribution is defined between 0 and 1.")
	}
	return math.Gamma(alpha+beta) / (math.Gamma(alpha) * math.Gamma(beta)) * math.Pow(x, alpha-1) * math.Pow(1-x, beta-1)
}

func GammaDist(x float64, alpha float64, beta float64) float64 {
	if alpha < 0 || beta < 0 || x < 0 {
		log.Fatalf("Alpha, beta parameters and input value must be greater than or equal to 0 in the gamma distribution.")
	}
	return (math.Pow(beta, alpha) / math.Gamma(alpha)) * math.Pow(x, alpha-1) * math.Exp(-beta*x)
}

func UninformativeGamma(x float64) float64 {
	return GammaDist(x, 1, 1)
}

//returns an instantiation of a normal distribution for a particular mean and SD
func NormalClosure(mu float64, sigma float64) func(float64) float64 {
	return func(x float64) float64 {
		return NormalDist(x, mu, sigma)
	}
}

func BetaClosure(alpha float64, beta float64) func(float64) float64 {
	return func(x float64) float64 {
		return BetaDist(x, alpha, beta)
	}
}

func GammaClosure(alpha float64, beta float64) func(float64) float64 {
	return func(x float64) float64 {
		return GammaDist(x, alpha, beta)
	}
}

func NormalLeftIntegral(x float64, mu float64, sigma float64) float64 {
	f := NormalClosure(mu, sigma)
	return DefiniteIntegral(f, mu-200*sigma, x)
}

func NormalRightIntegral(x float64, mu float64, sigma float64) float64 {
	f := NormalClosure(mu, sigma)
	return DefiniteIntegral(f, mu+200*sigma, x)
}

func NormalAdaptiveIntegral(left string, right string, mu float64, sigma float64) float64 {
	var leftInf, rightInf bool
	f := NormalClosure(mu, sigma)
	if left == "-INF" || left == "-Inf" || left == "-inf" {
		leftInf = true
	}
	if right == "INF" || right == "Inf" || right == "inf" {
		rightInf = true
	}
	if leftInf && rightInf {
		return 1.0
	} else if !leftInf && !rightInf {
		l := common.StringToFloat64(left)
		r := common.StringToFloat64(right)
		//if l > mu+10*sigma || r < mu-10*sigma {
		//	return DefiniteSmallIntegral(f, l, r)
		//}
		return DefiniteSmallIntegral(f, l, r)
	} else if leftInf {
		r := common.StringToFloat64(right)
		if r > mu+6*sigma { //Romberg can fail if a large right tail is evaluated in this case. R returns 1.0 for normal values over 6.
			return 1.0
		}
		if r < mu - 38*sigma {
			return 0.0
		}
		if r > mu-3 * sigma { //dealing with values close to mu
			return DefiniteSmallIntegral(f, r-15*sigma, r)
		} else {
			return DefiniteSmallIntegral(f, r-10*sigma, r)
		}
	} else if rightInf {
		l := common.StringToFloat64(left)
		if l < mu-6 * sigma {
			return 1.0 //same as above
		}
		if l > mu + 38 * sigma {
			return 0.0
		}
		if l < mu+10*sigma {
			return DefiniteSmallIntegral(f, l, l+15*sigma)
		} else {
			return DefiniteSmallIntegral(f, l, l+10*sigma)
		}
	} else {
		log.Fatalf("Something went wrong.")
		return -1.0
	}
}

func BetaLeftIntegral(x float64, alpha float64, beta float64) float64 {
	f := BetaClosure(alpha, beta)
	return DefiniteIntegral(f, 0, x)
}

func BetaRightIntegral(x float64, alpha float64, beta float64) float64 {
	f := BetaClosure(alpha, beta)
	return DefiniteIntegral(f, x, 1)
}

func BetaIntegral(left float64, right float64, alpha float64, beta float64) float64 {
	f := BetaClosure(alpha, beta)
	return DefiniteIntegral(f, left, right)
}

func GammaLeftIntegral(x float64, alpha float64, beta float64) float64 {
	f := GammaClosure(alpha, beta)
	return DefiniteIntegral(f, 0, x)
}

func GammaRightIntegral(x float64, alpha float64, beta float64) float64 {
	f := GammaClosure(alpha, beta)
	return 1 - DefiniteIntegral(f, 0, x)
}

func GammaIntegral(left float64, right float64, alpha float64, beta float64) float64 {
	f := GammaClosure(alpha, beta)
	return DefiniteIntegral(f, left, right)
}

func PoissonLeftSummation(k int, lambda float64) float64 {
	var answer float64 = 0
	for i := 0; i < k+1; i++ {
		answer = answer + PoissonDist(i, lambda)
	}
	return answer
}

//The Poisson right sum is infinite but can be evaluated as 1 - the left hand sum, which is finite.
func PoissonRightSummation(k int, lambda float64) float64 {
	return 1 - PoissonLeftSummation(k-1, lambda)
}

//inclusive on both ends
func PoissonSum(left int, right int, lambda float64) float64 {
	if right > left {
		log.Fatalf("PoissonSum failed. Right side value must be lower than the left side value.")
	}
	var answer float64 = 0
	for i := left; i < right; i++ {
		answer = answer + PoissonDist(i, lambda)
	}
	return answer
}

//all binomial summations include the query value n. Thus, the left and right sums will not add to 1.
func BinomialLeftSummation(n int, k int, p float64) float64 {
	if k <= n/2 {
		return evaluateLeftBinomialSum(n, k, p)
	} else {
		return 1 - evaluateRightBinomialSum(n, k+1, p)
	}
}

func BinomialRightSummation(n int, k int, p float64) float64 {
	if k > n/2 {
		return evaluateRightBinomialSum(n, k, p)
	} else {
		return 1 - evaluateLeftBinomialSum(n, k-1, p)
	}
}

//inclusive on both ends
func BinomialSum(left int, right int, n int, p float64) float64 {
	if right < left {
		log.Fatalf("BinomialSum failed. Right side value must be greater than the left side value.")
	}
	var answer float64 = 0
	for i := left; i <= right; i++ {
		answer = answer + BinomialDist(n, i, p)
	}
	return answer
}

func evaluateRightBinomialSum(n int, k int, p float64) float64 {
	var answer float64 = 0
	for i := k; i <= n; i++ {
		answer = answer + BinomialDist(n, i, p)
	}
	return answer
}

func evaluateLeftBinomialSum(n int, k int, p float64) float64 {
	var answer float64 = 0
	for i := 0; i < k+1; i++ {
		answer = answer + BinomialDist(n, i, p)
	}
	return answer
}

//TODO: Think of a way to test these two functions
//Measures the divergence between two probability distributions. Generally evaluated as an indefinite integral
//So set start and end to arbitrarily large numbers such that p(>end) -> 0 if the function is supported to infinity.
func ContinuousKullbackLeiblerDivergence(p func(float64) float64, q func(float64) float64, start float64, end float64) float64 {
	f := func(x float64) float64 {
		return p(x) * math.Log2(p(x)/q(x))
	}
	return DefiniteIntegral(f, start, end)
}

//inclusive range
func DiscreteKullbackLeiblerDivergence(p func(int) float64, q func(int) float64, start int, end int) float64 {
	f := func(x int) float64 {
		return p(x) * math.Log2(p(x)/q(x))
	}
	var answer float64 = 0
	for i := start; i < end+1; i++ {
		answer = answer + f(i)
	}
	return answer
}
