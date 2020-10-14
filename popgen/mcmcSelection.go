package popgen

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math"
	"math/rand"
	//DEBUG packages:
	"os"
	"runtime/pprof"
)

//To access debug prints, set verbose to 1 and then compile.
const verbose int = 1

type Theta struct {
	alpha       []float64
	mu          float64
	sigma       float64
	probability float64
	likelihood  float64
}

//MetropolisAccept is a helper function of MetropolisHastings that determines whether to accept or reject a candidate parameter set.
func MetropolisAccept(old Theta, thetaPrime Theta) bool {
	var pAccept, bayes, hastings, yRand float64
	yRand = math.Log(rand.Float64()) // random number from 0 to 1, transformed into log-space
	bayes = BayesRatio(old, thetaPrime)
	hastings = HastingsRatio(old, thetaPrime)
	pAccept = numbers.MultiplyLog(bayes, hastings)
	//pAccept = numbers.MinFloat64(1.0, BayesRatio(old, thetaPrime)*HastingsRatio(old, thetaPrime))
	if verbose > 0 {
		log.Printf("log(bayesRatio): %e, log(hastingsRatio): %e, log(likelihoodRatio): %e, log(rand): %e, decision: %t\n", bayes, hastings, pAccept, yRand, pAccept > yRand)
	}
	return pAccept > yRand
}

//HastingsRatio is a helper function of MetropolisAccept that returns the Hastings Ratio between two parameter sets.
func HastingsRatio(tOld Theta, tNew Theta) float64 {
	var newGivenOld, oldGivenNew float64
	newGivenOld = numbers.NormalDist(tNew.mu, tOld.mu, tOld.sigma) //* numbers.GammaDist(tNew.sigma, tOld.sigma*tOld.sigma, tOld.sigma)
	oldGivenNew = numbers.NormalDist(tOld.mu, tNew.mu, tNew.sigma) //* numbers.GammaDist(tOld.sigma, tNew.sigma*tNew.sigma, tNew.sigma)
	return math.Log(oldGivenNew / newGivenOld)
}

//BayesRatio is a helper function of MetropolisAccept taht returns the ratio of likelihoods of parameter sets
func BayesRatio(old Theta, thetaPrime Theta) float64 {
	like := numbers.DivideLog(thetaPrime.likelihood, old.likelihood)
	prob := math.Log(thetaPrime.probability / old.probability)
	if verbose > 0 {
		log.Printf("Old log(like): %e, New log(like): %e, likeRatio: %f, probRatio: %f\n", old.likelihood, thetaPrime.likelihood, math.Exp(like), math.Exp(prob))
	}
	return numbers.MultiplyLog(like, prob)
}

//GenerateCandidateThetaPrime is a helper function of Metropolis Hastings that picks a new set of parameters based on the state of the current parameter set t.
func GenerateCandidateThetaPrime(t Theta, data AFS, binomMap [][]float64) Theta {
	//sample from uninformative gamma
	var alphaPrime []float64
	var p float64 = 1.0
	var likelihood float64
	alphaPrime = make([]float64, len(t.alpha))

	//sample new sigma from a gamma function where the mean is always the current sigma value
	//mean of a gamma dist is alpha / beta, so mean = alpha / beta = sigma**2 / sigma = sigma
	//other condition is that the variance is fixed at 1 (var = alpha / beta**2 = sigma**2 / sigma**2
	//TODO: sigmaPrime still reverts to ultrasmall values, impeding step size. Need a permanant solution before this tool can be used effectively.
	//sigmaPrime := numbers.RandGamma(t.sigma*t.sigma, t.sigma)
	//sigmaPrime := numbers.RandGamma(1.0, 1.0)
	sigmaPrime := numbers.RandGamma(t.sigma*50, 50)
	//sigmaPrime = numbers.MaxFloat64(sigmaPrime, 0.01)
	//DEBUG: fmt.Printf("sigmaPrime: %e. tSigma: %e.\n", sigmaPrime, t.sigma)
	muPrime := numbers.SampleInverseNormal(t.mu, sigmaPrime)
	for i := 0; i < len(t.alpha); i++ {
		alphaPrime[i] = numbers.SampleInverseNormal(muPrime, sigmaPrime)
		p = p * numbers.NormalDist(alphaPrime[i], muPrime, sigmaPrime)
	}
	p = p * numbers.UninformativeGamma(sigmaPrime) * numbers.NormalDist(muPrime, t.mu, sigmaPrime)
	likelihood = AFSLikelihood(data, alphaPrime, binomMap)
	if verbose > 0 {
		log.Printf("Candidate Theta. Mu: %f. Sigma:%f. Probability:%e. Likelihood: %e.\n", muPrime, sigmaPrime, p, likelihood)
	}
	return Theta{alphaPrime, muPrime, sigmaPrime, p, likelihood}
}

//InitializeTheta is a helper function of Metropolis Hastings that generates the initial value of theta based on argument values.
func InitializeTheta(m float64, s float64, data AFS, binomMap [][]float64) Theta {
	k := len(data.sites)
	answer := Theta{mu: m, sigma: s}
	var p float64 = 1.0
	answer.alpha = make([]float64, k)
	for i := 0; i < k; i++ {
		answer.alpha[i] = numbers.SampleInverseNormal(m, s)
		p = p * numbers.NormalDist(answer.alpha[i], m, s)
	}
	//now multiply the probability of alpha, currently p, by the probability of drawing m and s from distributions if the previous state was m and s.
	answer.probability = p * numbers.UninformativeGamma(s) * numbers.NormalDist(m, m, s)
	answer.likelihood = AFSLikelihood(data, answer.alpha, binomMap)
	return answer
}

//MetropolisHastings implements the MH algorithm for Markov Chain Monte Carlo approximation of the posterior distribution for selection based on an input allele frequency spectrum.
//muZero and sigmaZero represent the starting hyperparameter values.
func MetropolisHastings(data AFS, muZero float64, sigmaZero float64, iterations int, outFile string) {
	//profiling test code DEBUG
	/*f, err := os.Create("testProfile.prof")
	if err != nil {
		log.Fatal(err)
	}
	pprof.StartCPUProfile(f)
	defer pprof.StopCPUProfile()*/

	out := fileio.EasyCreate(outFile)
	defer out.Close()

	if verbose > 0 {
		log.Println("Hello, I'm about to calculate MCMC.")
	}

	maxN := findMaxN(data)
	binomMap := make([][]float64, maxN+1)
	var currAccept bool
	if verbose > 0 {
		log.Println("Hello, I'm about to initialize theta.")
	}
	//initialization to uninformative standard normal
	t := InitializeTheta(muZero, sigmaZero, data, binomMap)
	if verbose > 0 {
		log.Printf("Initial Theta: mu: %f. sigma: %f. probability: %e. likelihood: %e.\n", t.mu, t.sigma, t.probability, t.likelihood)
	}
	fmt.Fprintf(out, "Iteration\tMu\tSigma\tAccept\n")

	for i := 0; i < iterations; i++ {
		tCandidate := GenerateCandidateThetaPrime(t, data, binomMap)
		if MetropolisAccept(t, tCandidate) {
			t = tCandidate
			currAccept = true
		} else {
			currAccept = false
		}
		fmt.Fprintf(out, "%v\t%e\t%e\t%t\n", i, t.mu, t.sigma, currAccept)
	}
}

//findMaxN is a helper function that aids in the generation of binomMap. In order to determine the proper length of the binomMap, we need to figure out which variant has the largest value of N.
func findMaxN(data AFS) int {
	var answer int = 0
	for i := 0; i < len(data.sites); i++ {
		answer = numbers.Max(answer, data.sites[i].n)
	}
	return answer
}
