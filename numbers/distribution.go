package numbers

import (
	"log"
	"math"
)

//returns the normal distribution value x for a distribution with mean mu and standard deviation sigma.
func NormalDist(x float64, mu float64, sigma float64) float64 {
	return (1 / (sigma * math.Sqrt(2*math.Pi))) * math.Exp(-0.5*math.Pow((x-mu)/sigma, 2))
}

func BinomialDist(n int, k int, p float64) float64 {
	return float64(BinomCoefficient(n, k)) * math.Pow(p, float64(k)) * math.Pow(1-p, float64(n-k))
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
