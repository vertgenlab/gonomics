package popgen

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/numbers"
	"math/rand"
	//DEBUG: "fmt"
)

type Theta struct {
	alpha       []float64
	mu          float64
	sigma       float64
	probability float64
}

func Metropolis_Accept(old Theta, thetaPrime Theta, data AFS, logSpace bool) bool {
	yRand := rand.Float64()
	var pAccept float64
	pAccept = common.MinFloat64(1.0, Bayes_Ratio(old, thetaPrime, data, logSpace)*Hastings_Ratio(old, thetaPrime))
	//DEBUG: fmt.Printf("Likelihood ratio: %f\n", pAccept)
	return pAccept > yRand
}

func Hastings_Ratio(tOld Theta, tNew Theta) float64 {
	var newGivenOld, oldGivenNew float64

	newGivenOld = numbers.NormalDist(tNew.mu, tOld.mu, tOld.sigma) //* numbers.GammaDist(tNew.sigma, tOld.sigma*tOld.sigma, tOld.sigma)
	oldGivenNew = numbers.NormalDist(tOld.mu, tNew.mu, tNew.sigma) //* numbers.GammaDist(tOld.sigma, tNew.sigma*tNew.sigma, tNew.sigma)

	/*
		for i := 0; i < len(tOld.alpha); i++ {
			newGivenOld = newGivenOld * numbers.NormalDist(tNew.alpha[i], tNew.mu, tNew.sigma)
			oldGivenNew = oldGivenNew * numbers.NormalDist(tOld.alpha[i], tOld.mu, tOld.sigma)
		}*/
	return oldGivenNew / newGivenOld
}

func Bayes_Ratio(old Theta, thetaPrime Theta, data AFS, logSpace bool) float64 {
	if logSpace {
		return AFSLogLikelihood(data, old.alpha) * old.probability / (AFSLogLikelihood(data, thetaPrime.alpha) * thetaPrime.probability)
	}
	return AFSLikelihood(data, old.alpha) * old.probability / (AFSLikelihood(data, thetaPrime.alpha) * thetaPrime.probability)
}

func GenerateCandidateThetaPrime(t Theta) Theta {
	//sample from uninformative gamma
	var alphaPrime []float64
	var p float64 = 1.0
	alphaPrime = make([]float64, len(t.alpha))

	//sample new sigma from a gamma function where the mean is always the current sigma value
	//mean of a gamma dist is alpha / beta, so mean = alpha / beta = sigma**2 / sigma = sigma
	//other condition is that the variance is fixed at 1 (var = alpha / beta**2 = sigma**2 / sigma**2
	//TODO: sigmaPrime still reverts to ultrasmall values, impeding step size. Need a permanant solution before this tool can be used effectively.
	//sigmaPrime := numbers.RandGamma(t.sigma*t.sigma, t.sigma)
	sigmaPrime := numbers.RandGamma(1.0, 1.0) 
	//sigmaPrime = common.MaxFloat64(sigmaPrime, 0.01)
	//DEBUG: fmt.Printf("sigmaPrime: %e. tSigma: %e.\n", sigmaPrime, t.sigma)
	muPrime := numbers.SampleInverseNormal(t.mu, sigmaPrime)
	for i := 0; i < len(t.alpha); i++ {
		alphaPrime[i] = numbers.SampleInverseNormal(muPrime, sigmaPrime)
		p = p * numbers.NormalDist(alphaPrime[i], muPrime, sigmaPrime)
	}
	p = p * numbers.UninformativeGamma(sigmaPrime) * numbers.NormalDist(muPrime, t.mu, sigmaPrime)

	return Theta{alphaPrime, muPrime, sigmaPrime, p}
}

func InitializeTheta(m float64, s float64, k int) Theta {
	answer := Theta{mu: m, sigma: s}
	var p float64 = 1.0
	answer.alpha = make([]float64, k)
	for i := 0; i < k; i++ {
		answer.alpha[i] = numbers.SampleInverseNormal(m, s)
		p = p * numbers.NormalDist(answer.alpha[i], m, s)
	}
	//now multiply the probability of alpha, currently p, by the probability of drawing m and s from distributions if the previous state was m and s.
	answer.probability = p * numbers.UninformativeGamma(s) * numbers.NormalDist(m, m, s)
	return answer
}

//MetropolisHastings implements the MH algorithm for Markov Chain Monte Carlo approximation of the posterior distribution for selection based on an input allele frequency spectrum.
//muZero and sigmaZero represent the starting hyperparameter values.
func MetropolisHastings(data AFS, muZero float64, sigmaZero float64, iterations int, logSpace bool) ([]float64, []float64, []bool) {
	muList := make([]float64, iterations)
	sigmaList := make([]float64, iterations)
	acceptList := make([]bool, iterations)
	//initialization to uninformative standard normal
	t := InitializeTheta(muZero, sigmaZero, len(data.sites))
	for i := 0; i < iterations; i++ {
		tCandidate := GenerateCandidateThetaPrime(t)
		if Metropolis_Accept(tCandidate, t, data, logSpace) {
			t = tCandidate
			acceptList[i] = true
		} else {
			acceptList[i] = false
		}
		muList[i] = t.mu
		sigmaList[i] = t.sigma
	}
	return muList, sigmaList, acceptList
}
