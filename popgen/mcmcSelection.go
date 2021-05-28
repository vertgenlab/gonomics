package popgen

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime/pprof"
)

//To access debug prints, set verbose to 1 or 2 and then compile. 2 returns lots of debug info, and 1 returns formatted debug info in tsv format for plotting.
const verbose int = 0

//The McmcSettings type stores settings for the various Mcmc helper functions.`
type McmcSettings struct {
	Iterations    int
	MuStep        float64
	MuZero        float64
	SigmaStep     float64
	SigmaZero     float64
	RandSeed      bool
	SetSeed       int64
	UnPolarized   bool
	DivergenceAscertainment bool
	FixedSigma	bool
	D	int //D is the size of the ascertainment subset.
	IntegralError float64
	Verbose int
}

//The Theta struct stores parameter sets, including the alpha vector, mu, and sigma parameters, along with the likelihood of a particular parameter set for MCMC.
type Theta struct {
	alpha []float64
	mu    float64
	sigma float64
	//probability float64
	likelihood float64
}

//MetropolisAccept is a helper function of MetropolisHastings that determines whether to accept or reject a candidate parameter set.
func MetropolisAccept(old Theta, thetaPrime Theta, s McmcSettings) bool {
	var pAccept, yRand float64
	yRand = math.Log(rand.Float64())
	var decision bool
	if thetaPrime.sigma < 0 || thetaPrime.sigma > 0.5 {//if sigma dips below zero or above 0.5, the candidate set is automatically discarded.
		return false
	}
	pAccept = BayesRatio(old, thetaPrime, s)
	decision = pAccept > yRand

	if s.Verbose == 1 {
		log.Printf("%e\t%e\t%e\t%e\t%e\t%e\t%t\n", old.mu, thetaPrime.mu, old.likelihood, thetaPrime.likelihood, pAccept, yRand, decision)
	}
	return decision
}
/*
//HastingsRatio is a helper function of MetropolisAccept that returns the Hastings Ratio (logspace) between two parameter sets.
func HastingsRatio(tOld Theta, tNew Theta, sigmaStep float64) float64 {
	var newGivenOld, oldGivenNew float64
	newGivenOld = numbers.GammaDist(tNew.sigma, sigmaStep, sigmaStep/tOld.sigma)
	oldGivenNew = numbers.GammaDist(tOld.sigma, sigmaStep, sigmaStep/tNew.sigma)
	return math.Log(oldGivenNew / newGivenOld)
}*/

//BayesRatio is a helper function of MetropolisAccept taht returns the ratio of likelihoods of parameter sets
func BayesRatio(old Theta, thetaPrime Theta, s McmcSettings) float64 {
	like := numbers.DivideLog(thetaPrime.likelihood, old.likelihood)
	//prob := numbers.DivideLog(thetaPrime.probability, old.probability)
	//prob = 0 //trick for debug
	if s.Verbose > 1 {
		log.Printf("Old log(like): %e, New log(like): %e, likeRatio: %f", old.likelihood, thetaPrime.likelihood, math.Exp(like))
	}
	return like
}

//GenerateCandidateThetaPrime is a helper function of Metropolis Hastings that picks a new set of parameters based on the state of the current parameter set t.
func GenerateCandidateThetaPrime(t Theta, data Afs, binomCache [][]float64, s McmcSettings) Theta {
	//sample from uninformative gamma
	var alphaPrime []float64
	//var p float64 = 0.0
	var likelihood, muPrime, sigmaPrime float64
	alphaPrime = make([]float64, len(t.alpha))
	sigmaPrime = numbers.SampleInverseNormal(t.sigma, s.SigmaStep)
	muPrime = numbers.SampleInverseNormal(t.mu, s.MuStep)
	for i := range t.alpha {
		alphaPrime[i] = numbers.SampleInverseNormal(muPrime, sigmaPrime)
	}

	if s.DivergenceAscertainment {
		likelihood = AfsDivergenceAscertainmentLikelihood(data, alphaPrime, binomCache, s.D, s.IntegralError)
	} else {
		likelihood = AfsLikelihood(data, alphaPrime, binomCache, s.IntegralError)
	}
	if s.Verbose > 1 {
		log.Printf("Candidate Theta. Mu: %f. Sigma:%f. LogLikelihood: %e.\n", muPrime, sigmaPrime, likelihood)
	}
	return Theta{alphaPrime, muPrime, sigmaPrime, likelihood}
}

//InitializeTheta is a helper function of Metropolis Hastings that generates the initial value of theta based on argument values.
func InitializeTheta(m float64, sig float64, data Afs, binomCache [][]float64, s McmcSettings) Theta {
	answer := Theta{mu: m, sigma: sig}
	answer.alpha = make([]float64, len(data.Sites))
	for i := range data.Sites {
		answer.alpha[i] = numbers.SampleInverseNormal(m, sig)
	}
	if s.DivergenceAscertainment {
		answer.likelihood = AfsDivergenceAscertainmentLikelihood(data, answer.alpha, binomCache, s.D, s.IntegralError)
	} else {
		answer.likelihood = AfsLikelihood(data, answer.alpha, binomCache, s.IntegralError)
	}
	return answer
}

//MetropolisHastings implements the MH algorithm for Markov Chain Monte Carlo approximation of the posterior distribution for selection based on an input allele frequency spectrum.
//muZero and sigmaZero represent the starting hyperparameter values.
func MetropolisHastings(data Afs, outFile string, s McmcSettings) {
	var err error
	if s.Verbose > 1 {
		f, err := os.Create("testProfile.prof")
		if err != nil {
			exception.PanicOnErr(err)
		}
		err = pprof.StartCPUProfile(f)
		if err != nil {
			exception.PanicOnErr(err)
		}
		defer pprof.StopCPUProfile()
	}

	out := fileio.EasyCreate(outFile)
	defer out.Close()

	if s.Verbose > 1 {
		log.Println("Hello, I'm about to calculate MCMC.")
	}
	allN := findAllN(data)
	binomCache := BuildBinomCache(allN)

	var currAccept bool
	if s.Verbose > 1 {
		log.Println("Hello, I'm about to initialize theta.")
	}
	//initialization to uninformative standard normal
	t := InitializeTheta(s.MuZero, s.SigmaZero, data, binomCache, s)
	if s.Verbose > 1 {
		log.Printf("Initial Theta: mu: %f. sigma: %f. LogLikelihood: %e.", t.mu, t.sigma, t.likelihood)
	}
	if s.Verbose == 1 {
		log.Printf("OldMu\tNewMu\tOldLikelihood\tNewLikelihood\tpAccept\tlogRand\tDecision\n")
	}
	_, err = fmt.Fprintf(out, "Iteration\tMu\tSigma\tAccept\n")
	if err != nil {
		exception.PanicOnErr(err)
	}
	for i := 0; i < s.Iterations; i++ {
		tCandidate := GenerateCandidateThetaPrime(t, data, binomCache, s)
		if MetropolisAccept(t, tCandidate, s) {
			t = tCandidate
			currAccept = true
		} else {
			currAccept = false
		}
		_, err = fmt.Fprintf(out, "%v\t%e\t%e\t%t\n", i, t.mu, t.sigma, currAccept)
		if err != nil {
			exception.PanicOnErr(err)
		}
	}
}

//BuildBinomCache makes a 2D matrix where each entry binomCache[n][k] is equal to [n choose k] in logSpace.
func BuildBinomCache(allN []int) [][]float64 {
	binomCache := make([][]float64, numbers.MaxIntSlice(allN)+1)

	var n, k int
	for n = range allN {
		binomCache[allN[n]] = make([]float64, allN[n])
		for k = 1; k < allN[n]; k++ {
			binomCache[allN[n]][k] = numbers.BinomCoefficientLog(allN[n], k)
		}
	}
	return binomCache
}

//findAllN is a helper function of Metropolis Hastings that returns all the unique values of N present in an input Afs struct.
func findAllN(data Afs) []int {
	var answer []int = make([]int, 0)
	for i := 0; i < len(data.Sites); i++ {
		if !common.IntSliceContains(answer, data.Sites[i].N) {
			answer = append(answer, data.Sites[i].N)
		}
	}
	return answer
}
