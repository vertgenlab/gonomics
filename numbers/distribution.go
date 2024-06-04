package numbers

import (
	"fmt"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"math"
)

// NormalDist returns the normal distribution value x for a distribution with mean mu and standard deviation sigma.
func NormalDist(x float64, mu float64, sigma float64) float64 {
	return (1 / (sigma * math.Sqrt(2*math.Pi))) * math.Exp(-0.5*math.Pow((x-mu)/sigma, 2))
}

// StandardNormalDist returns the probability density for an input x value on a standard normal distribution (mu=0, sigma=1).
func StandardNormalDist(x float64) float64 {
	return NormalDist(x, 0, 1)
}

// BinomialDist returns the probability mass from a binomial distribution with k successes out of n observations with success probability p.
// The second return is false if no overflow/underflow was detected. If underflow was detected, the program returns 0 and true.
// If logOutput is true, answer will be returned as log(answer).
func BinomialDist(n int, k int, p float64, logOutput bool) (float64, bool) {
	logAnswer := BinomialDistLog(n, k, p)
	if logOutput {
		return logAnswer, false
	}
	if logspace.CanConvert(logAnswer) {
		return math.Exp(logAnswer), false
	}
	return 0, true
}

// NegativeBinomialDist returns the probability mass of a negative binomial distribution with
// shape parameter r and success probability p at k. Pr(X = k)
// the second return is true if overflow/underflow was detected.
func NegativeBinomialDist(k int, r float64, p float64, logOutput bool) (float64, bool) {
	numerator, _ := math.Lgamma(float64(k) + r)
	denominatorLeft, _ := math.Lgamma(float64(k) + 1)
	denominatorRight, _ := math.Lgamma(r)
	coefficient := logspace.Divide(numerator, logspace.Multiply(denominatorLeft, denominatorRight))
	f := logspace.Pow(math.Log(1.0-p), float64(k))
	s := logspace.Pow(math.Log(p), r)
	var answer = logspace.Multiply(coefficient, logspace.Multiply(f, s))
	if logOutput {
		return answer, false
	}
	if !logspace.CanConvert(answer) {
		return 0, true
	}
	return math.Exp(answer), false
}

// GeometricDist returns the density of the geometric for k failures with success probability p.
// Note that this is the version of the geometric distribution with support from 0 to +INF.
func GeometricDist(k int, p float64) float64 {
	return math.Pow(1.0-p, float64(k)) * p
}

// ExpDist returns the density of the standard exponential distribution y=e^-x.
func ExpDist(x float64) float64 {
	return math.Exp(-x)
}

// PoissonDist returns the probability density of a poisson distribution with parameter lambda at the integer value k.
func PoissonDist(k int, lambda float64) float64 {
	if k < 0 {
		log.Fatalf("The poisson distribution is supported for k > 0.")
	}
	if lambda <= 0 {
		log.Fatalf("The poisson distribution is supported for lambda >= 0.")
	}
	return (math.Pow(lambda, float64(k)) * math.Pow(math.E, -lambda)) / float64(Factorial(k))
}

// BetaDist returns the probability density of a beta distribution with parameters alpha and beta at position x.
func BetaDist(x float64, alpha float64, beta float64) float64 {
	if alpha <= 0 {
		log.Fatalf("Alpha parameter must be greater than 0.")
	}
	if beta <= 0 {
		log.Fatalf("Beta parameter must be greater than 0. Beta: %f.", beta)
	}
	if x < 0 || x > 1 {
		log.Fatalf("Value x out of range. The beta distribution is defined between 0 and 1.")
	}
	return math.Pow(x, alpha-1) * math.Pow(1-x, beta-1) / BetaFunc(alpha, beta)
}

// BetaFunc returns B(x, y), where B is the Beta Function, also known as the Euler integral of the first kind.
func BetaFunc(x float64, y float64) float64 {
	return math.Gamma(x) * math.Gamma(y) / math.Gamma(x+y)
}

// GammaDist returns the probability density of a gamma distribution with parameters alpha and beta at position x.
// alpha is the shape parameter and beta is the rate parameter.
func GammaDist(x float64, alpha float64, beta float64) float64 {
	if alpha < 0 || beta < 0 || x < 0 {
		log.Fatalf("Alpha, beta parameters and input value must be greater than or equal to 0 in the gamma distribution.")
	}
	return (math.Pow(beta, alpha) / math.Gamma(alpha)) * math.Pow(x, alpha-1) * math.Exp(-beta*x)
}

// NormalClosure returns an instantiation of a normal distribution for a particular mean mu and standard deviation sigma.
func NormalClosure(mu float64, sigma float64) func(float64) float64 {
	return func(x float64) float64 {
		return NormalDist(x, mu, sigma)
	}
}

// BetaClosure returns an instantiation of a Beta Distribution with fixed alpha and beta parameters.
func BetaClosure(alpha float64, beta float64) func(float64) float64 {
	return func(x float64) float64 {
		return BetaDist(x, alpha, beta)
	}
}

// GamaClosure returns an instantiation of a Gamma Distribution with fixed alpha and beta parameters.
func GammaClosure(alpha float64, beta float64) func(float64) float64 {
	return func(x float64) float64 {
		return GammaDist(x, alpha, beta)
	}
}

// NormalLeftIntegral returns the area under the curve of an input normal probability distribution defined by mean (mu) and standard deviation (sigma) to the left of an input point x.
func NormalLeftIntegral(x float64, mu float64, sigma float64) float64 {
	f := NormalClosure(mu, sigma)
	return DefiniteIntegral(f, mu-200*sigma, x)
}

// NormalRightIntegral returns the area under the curve of an input normal probability distribution defined by mean (mu) and standard deviation (sigma) to the right of an input point x.
func NormalRightIntegral(x float64, mu float64, sigma float64) float64 {
	f := NormalClosure(mu, sigma)
	return DefiniteIntegral(f, mu+200*sigma, x)
}

func LogNormalRightTailCDF(x, mu, sigma float64) (float64, error) {
	z := (x - mu) / sigma
	logErfc := math.Log(math.Erfc(z / math.Sqrt2))
	logHalf := math.Log(0.5)
	result := logHalf + logErfc

	if math.IsInf(result, 0) {
		return result, fmt.Errorf("overflow detected")
	}
	if result == math.Inf(-1) {
		return result, fmt.Errorf("underflow detected")
	}
	return result, nil
}

// NormalAdaptiveIntegral returns the integral under a normal probability distribution with mean mu and standard deviation sigma from a specified left and right bound.
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
		l := parse.StringToFloat64(left)
		r := parse.StringToFloat64(right)
		//if l > mu+10*sigma || r < mu-10*sigma {
		//	return DefiniteSmallIntegral(f, l, r)
		//}
		return DefiniteSmallIntegral(f, l, r)
	} else if leftInf {
		r := parse.StringToFloat64(right)
		if r > mu+6*sigma { //Romberg can fail if a large right tail is evaluated in this case. R returns 1.0 for normal values over 6.
			return 1.0
		}
		if r < mu-38*sigma {
			return 0.0
		}
		if r > mu-3*sigma { //dealing with values close to mu
			return DefiniteSmallIntegral(f, r-15*sigma, r)
		} else {
			return DefiniteSmallIntegral(f, r-10*sigma, r)
		}
	} else if rightInf {
		l := parse.StringToFloat64(left)
		if l < mu-6*sigma {
			return 1.0 //same as above
		}
		if l > mu+38*sigma {
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

// BetaLeftIntegral returns the integral to the left of an input point x from a beta distribution with parameters alpha and beta.
func BetaLeftIntegral(x float64, alpha float64, beta float64) float64 {
	f := BetaClosure(alpha, beta)
	return DefiniteIntegral(f, 0, x)
}

// BetaRightIntegral returns the integral to the right of an input point x from a beta distribution with parameters alpha and beta.
func BetaRightIntegral(x float64, alpha float64, beta float64) float64 {
	f := BetaClosure(alpha, beta)
	return DefiniteIntegral(f, x, 1)
}

// BetaIntegral calculates the integral under a beta distribution with parameters alpha and beta between a specified left and right bound.
func BetaIntegral(left float64, right float64, alpha float64, beta float64) float64 {
	f := BetaClosure(alpha, beta)
	return DefiniteIntegral(f, left, right)
}

// GammaLeftIntegral calculates the integral to the left of an input position on a gamma distribution with parameters alpha and beta.
func GammaLeftIntegral(x float64, alpha float64, beta float64) float64 {
	f := GammaClosure(alpha, beta)
	return DefiniteIntegral(f, 0, x)
}

// GammaRightIntegral calculates the integral to the right of an input position on a gamma distribution with parameters alpha and beta.
func GammaRightIntegral(x float64, alpha float64, beta float64) float64 {
	f := GammaClosure(alpha, beta)
	return 1 - DefiniteIntegral(f, 0, x)
}

// GammaIntegral calculates the integral between an input left and right bound of a gamma distribution with parameters alpha and beta.
func GammaIntegral(left float64, right float64, alpha float64, beta float64) float64 {
	f := GammaClosure(alpha, beta)
	return DefiniteIntegral(f, left, right)
}

// PoissonLeftSummation calculates the sum of probabilities to the left of an integer k, inclusive.
func PoissonLeftSummation(k int, lambda float64) float64 {
	var answer float64 = 0
	for i := 0; i < k+1; i++ {
		answer = answer + PoissonDist(i, lambda)
	}
	return answer
}

// PoissonRightSummation calculates the sum of probabilities to the right of an integer k, inclusive.
func PoissonRightSummation(k int, lambda float64) float64 {
	return 1 - PoissonLeftSummation(k-1, lambda)
}

// PoissonSum calculates the sum of probabilities between an input left and right bound, inclusive on both ends.
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

// BinomialLeftSummation calculates the sum of binomial probabilities to the left of k successes for a binomial distribution with n experiments and a success probability of p, inclusive.
func BinomialLeftSummation(n int, k int, p float64, logOutput bool) float64 {
	if n == k {
		if logOutput {
			return 0
		} else {
			return 1
		}
	}
	return evaluateLeftBinomialSum(n, k, p, logOutput)
}

// BinomialRightSummation calculates the sum of binomial probabilities to the right of k successes for a binomial distribution with n experiments and a success probability of p, inclusive.
func BinomialRightSummation(n int, k int, p float64, logOutput bool, exact bool) float64 {
	if k == 0 {
		if logOutput {
			return 0
		} else {
			return 1
		}
	}
	if float64(n)*p > 5 && float64(n)*(1-p) > 5 && !exact {
		return evaluateRightBinomialSumApproximate(n, k, p, logOutput)
	} else {
		log.Print("ended up here")
		return evaluateRightBinomialSum(n, k, p, logOutput)
	}
}

// BinomialSum calculates the sum of probabilities in a binomial distribution with n experiments and success probability p between two input k values. Inclusive on both ends.
func BinomialSum(left int, right int, n int, p float64, logOutput bool) float64 {
	if right < left {
		log.Fatalf("BinomialSum failed. Right side value must be greater than the left side value.")
	}
	var answer, _ = BinomialDist(n, left, p, logOutput)
	var curr float64
	for i := left; i <= right; i++ {
		curr, _ = BinomialDist(n, i, p, logOutput)
		if logOutput {
			answer = logspace.Add(answer, curr)
		} else {
			answer += curr
		}
	}
	return answer
}

// evaluateRightBinomialSum is a helper function that calculates the sum of probabilities under a binomial distribution with parameters n and p to the right of an input k value, inclusive.
func evaluateRightBinomialSum(n int, k int, p float64, logOutput bool) float64 {
	var curr float64
	var answer, _ = BinomialDist(n, k, p, logOutput)
	for i := k + 1; i <= n; i++ {
		curr, _ = BinomialDist(n, i, p, logOutput)
		if logOutput {
			answer = logspace.Add(answer, curr)
		} else {
			answer += curr
		}
	}
	return answer
}

// evaluateLeftBinomialSum is a helper function that calculates the sum of probabilities under a binomial distribution with parameters n and p to the left of an input k value, inclusive.
func evaluateLeftBinomialSum(n int, k int, p float64, logOutput bool) float64 {
	var curr float64
	var answer, _ = BinomialDist(n, k, p, logOutput)
	for i := 0; i < k; i++ {
		curr, _ = BinomialDist(n, i, p, logOutput)
		if logOutput {
			answer = logspace.Add(answer, curr)
		} else {
			answer += curr
		}
	}
	return answer
}

// evaluateRightBinomialSumApproximate is a helper function that calculates the sum of probabilities under a binomial distribution with parameters n and p to the right of an input k value, inclusive.
func evaluateRightBinomialSumApproximate(n int, k int, p float64, logOutput bool) float64 {
	var curr float64
	var x, mu, sig float64
	mu = float64(n) * p
	x = float64(k) - 0.5
	sig = math.Sqrt(float64(n) * p * (1 - p))

	var answer float64
	if logOutput {
		answer, _ = LogNormalRightTailCDF(x, mu, sig)
	} else {
		answer = NormalDist(x, mu, sig)
		for i := int(x) + 1; i <= n; i++ {
			curr = NormalDist(float64(i), mu, sig)
			answer += curr
		}
	}

	return answer
}
