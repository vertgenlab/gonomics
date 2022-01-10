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
	Iterations              int
	MuStep                  float64
	MuZero                  float64
	SigmaStep               float64
	SigmaZero               float64
	RandSeed                bool  //pseudorandom number generation control: set the seed as current unix time
	SetSeed                 int64 //pseudorandom number generation control: set the seed to a user-specified int64.
	UnPolarized             bool
	DivergenceAscertainment bool
	FixedSigma              bool
	D                       int //D is the size of the ascertainment subset.
	IntegralError           float64
	Verbose                 int
}

// The Theta struct stores parameter sets, including the alpha vector, mu, and sigma parameters, along with the likelihood of a particular parameter set for MCMC.
type Theta struct {
	alpha      []float64
	mu         float64
	sigma      float64
	prior      float64
	likelihood float64
}

// MetropolisAccept is a helper function of MetropolisHastings that determines whether to accept or reject a candidate parameter set.
func MetropolisAccept(old Theta, thetaPrime Theta, s McmcSettings) bool {
	var pAccept, yRand float64
	yRand = math.Log(rand.Float64())
	var decision bool

	/* I think the prior calculations can now handle this
	if thetaPrime.sigma < 0 || thetaPrime.sigma > 0.5 { //if sigma dips below zero or above 0.5, the candidate set is automatically discarded.
		return false
	}*/
	pAccept = PosteriorOdds(old, thetaPrime)
	decision = pAccept > yRand

	if s.Verbose == 1 {
		log.Printf("%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%t\n", old.mu, thetaPrime.mu, old.sigma, thetaPrime.sigma, old.likelihood, thetaPrime.likelihood, old.prior, thetaPrime.prior, pAccept, yRand, decision)
		//log.Printf("%g\t%g\t%g\t%g\t%g\t%g\t%t\n", old.mu, thetaPrime.mu, math.Exp(numbers.DivideLog(old.likelihood, thetaPrime.likelihood)), math.Exp(numbers.DivideLog(old.prior, thetaPrime.prior)), math.Exp(pAccept), math.Exp(yRand), decision)
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

// PosteriorOdds is a helper function of MetropolisAccept that returns the Bayes Factor times the Prior Odds
// this should be the probability of accepting (can be greater than 1) if the Hastings Ratio is one.
func PosteriorOdds(old Theta, thetaPrime Theta) float64 {
	bayesFactor := numbers.DivideLog(thetaPrime.likelihood, old.likelihood)
	priorOdds := numbers.DivideLog(thetaPrime.prior, old.prior)
	posteriorOdds := numbers.MultiplyLog(bayesFactor, priorOdds)
	return posteriorOdds
}

// prior returns log(probability) of having meanAlpha and sigma as mean
// and standard deviation of the function that will be generating the individual
// alpha values
func priorProb(mu float64, sigma float64) float64 {
	var sigmaPrior, muPrior float64

	//prior on sigma is a uniform distribution between zero and 0.5
	if sigma < 0 {
		return math.Inf(-1) // prior probability is zero
	} else {
		sigmaPrior = numbers.GammaDist(sigma, 2, 10)
	}

	//prior on alpha is normal with a mean of zero and a stdev of 3
	muPrior = numbers.NormalDist(mu, 0, 3)

	return math.Log(muPrior * sigmaPrior)
}

// GenerateCandidateThetaPrime is a helper function of Metropolis Hastings that picks a new set of
// parameters based on the state of the current parameter set t.
// TODO: We could avoid some memory allocations by passing in an "old" theta and overwriting the values
func GenerateCandidateThetaPrime(t Theta, data Afs, binomCache [][]float64, s McmcSettings) Theta {
	//sample from uninformative gamma
	var alphaPrime []float64
	var prior, likelihood, muPrime, sigmaPrime float64
	alphaPrime = make([]float64, len(t.alpha))

	if s.FixedSigma {
		sigmaPrime = t.sigma
	} else {
		sigmaPrime = numbers.SampleInverseNormal(t.sigma, s.SigmaStep)
	}

	muPrime = numbers.SampleInverseNormal(t.mu, s.MuStep)
	for i := range t.alpha {
		alphaPrime[i] = numbers.SampleInverseNormal(muPrime, sigmaPrime)
	}

	if s.DivergenceAscertainment {
		likelihood = AfsDivergenceAscertainmentLikelihood(data, alphaPrime, binomCache, s.D, s.IntegralError)
	} else {
		likelihood = AfsLikelihood(data, alphaPrime, binomCache, s.IntegralError)
	}

	prior = priorProb(muPrime, sigmaPrime)

	if s.Verbose > 1 {
		log.Printf("Candidate Theta. Mu: %f. Sigma:%f. LogLikelihood: %e.\n", muPrime, sigmaPrime, likelihood)
	}
	return Theta{alphaPrime, muPrime, sigmaPrime, prior, likelihood}
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
	answer.prior = priorProb(answer.mu, answer.sigma)
	return answer
}

//MetropolisHastings implements the MH algorithm for Markov Chain Monte Carlo approximation of the posterior distribution for selection based on an input allele frequency spectrum.
//muZero and sigmaZero represent the starting hyperparameter values.
func MetropolisHastings(data Afs, outFile string, s McmcSettings) {
	var err error
	if s.Verbose > 1 {
		f, err := os.Create("testProfile.prof")
		exception.PanicOnErr(err)
		err = pprof.StartCPUProfile(f)
		exception.PanicOnErr(err)
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
		log.Printf("OldMu\tNewMu\tOldLikelihood\tNewLikelihood\tOldPrior\tNewPrior\tpAccept\tlogRand\tDecision\n")
	}
	_, err = fmt.Fprintf(out, "Iteration\tMu\tSigma\tAccept\n")
	exception.PanicOnErr(err)
	for i := 0; i < s.Iterations; i++ {
		tCandidate := GenerateCandidateThetaPrime(t, data, binomCache, s)
		if MetropolisAccept(t, tCandidate, s) {
			t = tCandidate
			currAccept = true
		} else {
			currAccept = false
		}
		_, err = fmt.Fprintf(out, "%v\t%e\t%e\t%t\n", i, t.mu, t.sigma, currAccept)
		exception.PanicOnErr(err)
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
