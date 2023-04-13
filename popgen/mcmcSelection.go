package popgen

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime/pprof"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"golang.org/x/exp/slices"
)

// To access debug prints, set verbose to 1 or 2 and then compile. 2 returns lots of debug info, and 1 returns formatted debug info in tsv format for plotting.
const verbose int = 0

// The McmcSettings type stores settings for the various Mcmc helper functions.`
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
	SigmaPriorAlpha         float64 //defines the alpha parameter for the prior distribution for the Theta hyperparameter sigma.
	SigmaPriorBeta          float64 //defines the beta parameter for the prior distribution for the Theta hyperparameter sigma.
	MuPriorMean             float64 //defines the mean of the prior distribution for the Theta hyperparameter mu.
	MuPriorSigma            float64 //defines the standard deviation of the prior distribution for the Theta hyperparameter mu.
	IncludeRef              bool    //If true, includes the reference genome allele state as a datapoint in the allele frequency spectrum.
}

// The Theta struct stores parameter sets, including the alpha vector, mu, and sigma parameters, along with the likelihood of a particular parameter set for MCMC.
type Theta struct {
	alpha        []float64 //defines the vector of selction parameters alpha for each segregating site.
	mu           float64   //hyperparameter to generate alpha. Defines the mean of the distribution of alpha values.
	sigma        float64   //hyperparameter to generate alpha. Defines the standard deviation of the distribution of alpha values.
	priorDensity float64   //density of the prior distribution for a particular theta set.
	likelihood   float64   //likelihood value for a particular Theta set.
}

// MetropolisAccept is a helper function of MetropolisHastings that determines whether to accept or reject a candidate parameter set.
func MetropolisAccept(old Theta, thetaPrime Theta, s McmcSettings) bool {
	var pAccept, yRand float64
	yRand = math.Log(rand.Float64())
	var decision bool
	pAccept = PosteriorOdds(old, thetaPrime)
	decision = pAccept > yRand

	if s.Verbose == 1 {
		log.Printf("%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%t\n", old.mu, thetaPrime.mu, old.sigma, thetaPrime.sigma, old.likelihood, thetaPrime.likelihood, old.priorDensity, thetaPrime.priorDensity, pAccept, yRand, decision)
	}
	return decision
}

// PosteriorOdds is a helper function of MetropolisAccept that returns the Bayes Factor times the Prior Odds
// this should be the probability of accepting (can be greater than 1) if the Hastings Ratio is one.
func PosteriorOdds(old Theta, thetaPrime Theta) float64 {
	if thetaPrime.priorDensity == math.Inf(-1) || thetaPrime.likelihood == math.Inf(-1) { //avoid divide by -Inf error when the candidate set is overdispersed.
		return math.Inf(-1)
	}
	bayesFactor := logspace.Divide(thetaPrime.likelihood, old.likelihood)
	priorOdds := logspace.Divide(thetaPrime.priorDensity, old.priorDensity)
	posteriorOdds := logspace.Multiply(bayesFactor, priorOdds)
	return posteriorOdds
}

// priorProb returns log(probability) of having meanAlpha and sigma as mean
// and standard deviation of the function that will be generating the individual
// alpha values
func priorProb(mu float64, sigma float64, s McmcSettings) float64 {
	var sigmaPrior, muPrior float64

	if sigma < 0 {
		return math.Inf(-1) // prior probability is zero
	} else {
		sigmaPrior = numbers.GammaDist(sigma, s.SigmaPriorAlpha, s.SigmaPriorBeta)
	}
	muPrior = numbers.NormalDist(mu, s.MuPriorMean, s.MuPriorSigma)
	return math.Log(muPrior * sigmaPrior)
}

// GenerateCandidateThetaPrime is a helper function of Metropolis Hastings that picks a new set of
// parameters based on the state of the current parameter set t.
// TODO: We could avoid some memory allocations by passing in an "old" theta and overwriting the values
func GenerateCandidateThetaPrime(t Theta, data Afs, binomCache [][]float64, s McmcSettings) Theta {
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
	prior = priorProb(muPrime, sigmaPrime, s)

	if prior == math.Inf(-1) {
		likelihood = math.Inf(-1)
	} else if s.DivergenceAscertainment {
		likelihood = AfsDivergenceAscertainmentLikelihood(data, alphaPrime, binomCache, s.D, s.IntegralError)
	} else {
		likelihood = AfsLikelihood(data, alphaPrime, binomCache, s.IntegralError)
	}

	if s.Verbose > 1 {
		log.Printf("Candidate Theta. Mu: %f. Sigma:%f. LogLikelihood: %e.\n", muPrime, sigmaPrime, likelihood)
	}
	return Theta{alphaPrime, muPrime, sigmaPrime, prior, likelihood}
}

// InitializeTheta is a helper function of Metropolis Hastings that generates the initial value of theta based on argument values.
func InitializeTheta(m float64, sig float64, data Afs, binomCache [][]float64, s McmcSettings) Theta {
	answer := Theta{mu: m, sigma: sig}
	answer.alpha = make([]float64, len(data.Sites))
	for i := range data.Sites {
		answer.alpha[i] = numbers.SampleInverseNormal(m, sig)
	}
	answer.priorDensity = priorProb(answer.mu, answer.sigma, s)
	if answer.priorDensity == math.Inf(-1) {
		log.Fatalf("Initial theta set is too overdispersed to have a finite prior density in logSpace.")
	} else if s.DivergenceAscertainment {
		answer.likelihood = AfsDivergenceAscertainmentLikelihood(data, answer.alpha, binomCache, s.D, s.IntegralError)
	} else {
		answer.likelihood = AfsLikelihood(data, answer.alpha, binomCache, s.IntegralError)
	}

	return answer
}

// MetropolisHastings implements the MH algorithm for Markov Chain Monte Carlo approximation of the posterior distribution for selection based on an input allele frequency spectrum.
// muZero and sigmaZero represent the starting hyperparameter values.
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

	if s.Verbose > 1 {
		log.Println("Hello, I'm about to calculate MCMC.")
	}
	allN := findAllN(data)
	binomCache := BuildBinomCache(allN)

	var currAccept bool
	if s.Verbose > 1 {
		log.Println("Hello, I'm about to initialize theta.")
	}
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

	err = out.Close()
	exception.PanicOnErr(err)
}

// BuildBinomCache makes a 2D matrix where each entry binomCache[n][k] is equal to [n choose k] in logSpace.
func BuildBinomCache(allN []int) [][]float64 {
	binomCache := make([][]float64, numbers.MaxMany(allN...)+1)

	var n, k int
	for n = range allN {
		binomCache[allN[n]] = make([]float64, allN[n])
		for k = 1; k < allN[n]; k++ {
			binomCache[allN[n]][k] = numbers.BinomCoefficientLog(allN[n], k)
		}
	}
	return binomCache
}

// findAllN is a helper function of Metropolis Hastings that returns all the unique values of N present in an input Afs struct.
func findAllN(data Afs) []int {
	var answer = make([]int, 0)
	for i := 0; i < len(data.Sites); i++ {
		if !slices.Contains(answer, data.Sites[i].N) {
			answer = append(answer, data.Sites[i].N)
		}
	}
	return answer
}
