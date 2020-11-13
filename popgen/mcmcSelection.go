package popgen

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime/pprof"
)

//To access debug prints, set verbose to 1 or 2 and then compile. 2 returns lots of debug info, and 1 returns formatted debug info in tsv format for plotting.
const verbose int = 1
const step float64 = 50.0

type Theta struct {
	alpha []float64
	mu    float64
	sigma float64
	//probability float64
	likelihood float64
}

//MetropolisAccept is a helper function of MetropolisHastings that determines whether to accept or reject a candidate parameter set.
func MetropolisAccept(old Theta, thetaPrime Theta) bool {
	var pAccept, bayes, hastings, yRand float64
	yRand = math.Log(rand.Float64())
	var decision bool
	bayes = BayesRatio(old, thetaPrime)
	hastings = HastingsRatio(old, thetaPrime)
	pAccept = numbers.MultiplyLog(bayes, hastings)
	decision = pAccept > yRand
	//pAccept = numbers.MinFloat64(1.0, BayesRatio(old, thetaPrime)*HastingsRatio(old, thetaPrime))
	if verbose > 1 {
		log.Printf("bayesRatio: %e, hastingsRatio: %e, log(likelihoodRatio): %e, log(rand): %e, decision: %t\n", bayes, hastings, pAccept, yRand, decision)
	}
	if verbose == 1 {
		log.Printf("%e\t%e\t%e\t%e\t%e\t%e\t%t\n", old.mu, thetaPrime.mu, old.likelihood, thetaPrime.likelihood, pAccept, yRand, decision)
	}
	return decision
}

//HastingsRatio is a helper function of MetropolisAccept that returns the Hastings Ratio (logspace) between two parameter sets.
func HastingsRatio(tOld Theta, tNew Theta) float64 {
	var newGivenOld, oldGivenNew float64
	newGivenOld = numbers.NormalDist(tNew.mu, tOld.mu, tOld.sigma) * numbers.GammaDist(tNew.sigma, step, step/tOld.sigma)
	oldGivenNew = numbers.NormalDist(tOld.mu, tNew.mu, tNew.sigma) * numbers.GammaDist(tOld.sigma, step, step/tNew.sigma)
	return math.Log(oldGivenNew / newGivenOld)
}

//BayesRatio is a helper function of MetropolisAccept taht returns the ratio of likelihoods of parameter sets
func BayesRatio(old Theta, thetaPrime Theta) float64 {
	like := numbers.DivideLog(thetaPrime.likelihood, old.likelihood)
	//prob := numbers.DivideLog(thetaPrime.probability, old.probability)
	//prob = 0 //trick for debug
	if verbose > 1 {
		log.Printf("Old log(like): %e, New log(like): %e, likeRatio: %f", old.likelihood, thetaPrime.likelihood, math.Exp(like))
	}
	return like
}

//GenerateCandidateThetaPrime is a helper function of Metropolis Hastings that picks a new set of parameters based on the state of the current parameter set t.
func GenerateCandidateThetaPrime(t Theta, data AFS, binomCache [][]float64) Theta {
	//sample from uninformative gamma
	var alphaPrime []float64
	//var p float64 = 0.0
	var likelihood float64
	alphaPrime = make([]float64, len(t.alpha))

	//sample new sigma from a gamma function where the mean is always the current sigma value
	//mean of a gamma dist is alpha / beta, so mean = alpha / beta = sigma**2 / sigma = sigma
	//other condition is that the variance is fixed at 1 (var = alpha / beta**2 = sigma**2 / sigma**2
	sigmaPrime, _ := numbers.RandGamma(step, step/t.sigma)
	muPrime := numbers.SampleInverseNormal(t.mu, sigmaPrime)
	for i := 0; i < len(t.alpha); i++ {
		alphaPrime[i] = numbers.SampleInverseNormal(muPrime, sigmaPrime)
		//p = p * numbers.NormalDist(alphaPrime[i], muPrime, sigmaPrime)
		//p = numbers.MultiplyLog(p, math.Log(numbers.NormalDist(alphaPrime[i], muPrime, sigmaPrime)))
	}
	//p = numbers.MultiplyLog(p, math.Log(numbers.NormalDist(muPrime, t.mu, sigmaPrime)))
	//p = numbers.MultiplyLog(p, math.Log(numbers.UninformativeGamma(sigmaPrime)))
	likelihood = AFSLikelihood(data, alphaPrime, binomCache)
	if verbose > 1 {
		log.Printf("Candidate Theta. Mu: %f. Sigma:%f. LogLikelihood: %e.\n", muPrime, sigmaPrime, likelihood)
	}
	return Theta{alphaPrime, muPrime, sigmaPrime, likelihood}
}

//InitializeTheta is a helper function of Metropolis Hastings that generates the initial value of theta based on argument values.
func InitializeTheta(m float64, s float64, data AFS, binomCache [][]float64) Theta {
	k := len(data.sites)
	answer := Theta{mu: m, sigma: s}
	//var p float64 = 0.0
	answer.alpha = make([]float64, k)
	for i := 0; i < k; i++ {
		answer.alpha[i] = numbers.SampleInverseNormal(m, s)
		//p = p * numbers.NormalDist(answer.alpha[i], m, s)
		//	p = numbers.MultiplyLog(p, math.Log(numbers.NormalDist(answer.alpha[i], m, s)))
	}
	//now multiply the probability of alpha, currently p, by the probability of drawing m and s from distributions if the previous state was m and s.
	//answer.probability = p * numbers.UninformativeGamma(s) * numbers.NormalDist(m, m, s)
	//answer.probability = numbers.MultiplyLog(p,  math.Log(numbers.UninformativeGamma(s)))
	//answer.probability = numbers.MultiplyLog(p, math.Log(numbers.NormalDist(m, m, s)))
	answer.likelihood = AFSLikelihood(data, answer.alpha, binomCache)
	return answer
}

//MetropolisHastings implements the MH algorithm for Markov Chain Monte Carlo approximation of the posterior distribution for selection based on an input allele frequency spectrum.
//muZero and sigmaZero represent the starting hyperparameter values.
func MetropolisHastings(data AFS, muZero float64, sigmaZero float64, iterations int, outFile string) {
	if verbose > 1 {
		f, err := os.Create("testProfile.prof")
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	out := fileio.EasyCreate(outFile)
	defer out.Close()

	if verbose > 1 {
		log.Println("Hello, I'm about to calculate MCMC.")
	}

	maxN := findMaxN(data)
	binomCache := make([][]float64, maxN+1)

	allN := findAllN(data)
	var n, k int
	for n = 0; n < len(allN); n++ {
		binomCache[allN[n]] = make([]float64, allN[n])
		for k = 1; k < allN[n]; k++ {
			binomCache[allN[n]][k] = numbers.BinomCoefficientLog(allN[n], k)
		}
	}

	var currAccept bool
	if verbose > 1 {
		log.Println("Hello, I'm about to initialize theta.")
	}
	//initialization to uninformative standard normal
	t := InitializeTheta(muZero, sigmaZero, data, binomCache)
	if verbose > 1 {
		log.Printf("Initial Theta: mu: %f. sigma: %f. LogLikelihood: %e.", t.mu, t.sigma, t.likelihood)
	}
	if verbose == 1 {
		log.Printf("OldMu\tNewMu\tOldLikelihood\tNewLikelihood\tpAccept\tlogRand\tDecision\n")
	}
	fmt.Fprintf(out, "Iteration\tMu\tSigma\tAccept\n")

	for i := 0; i < iterations; i++ {
		tCandidate := GenerateCandidateThetaPrime(t, data, binomCache)
		if MetropolisAccept(t, tCandidate) {
			t = tCandidate
			currAccept = true
		} else {
			currAccept = false
		}
		fmt.Fprintf(out, "%v\t%e\t%e\t%t\n", i, t.mu, t.sigma, currAccept)
	}
}

//findMaxN is a helper function that aids in the generation of binomCache. In order to determine the proper length of the binomCache, we need to figure out which variant has the largest value of N.
func findMaxN(data AFS) int {
	var answer int = 0
	for i := 0; i < len(data.sites); i++ {
		answer = numbers.Max(answer, data.sites[i].n)
	}
	return answer
}

func findAllN(data AFS) []int {
	var answer []int = make([]int, 0)
	for i := 0; i < len(data.sites); i++ {
		if !common.IntSliceContains(answer, data.sites[i].n) {
			answer = append(answer, data.sites[i].n)
		}
	}
	return answer
}
