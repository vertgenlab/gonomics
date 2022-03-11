// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/popgen"
	"log"
	"math/rand"
)

func selectionMcmc(filename string, outFile string, s popgen.McmcSettings) {
	rand.Seed(s.SetSeed)
	data, err := popgen.VcfToAfs(filename, s.UnPolarized, s.DivergenceAscertainment, s.IncludeRef) //VcfToAFS is written with polarized as the argument for clarity, so the bool is flipped here.
	exception.FatalOnErr(err)
	popgen.MetropolisHastings(*data, outFile, s)
}

func usage() {
	fmt.Print(
		"selectionMCMC - Returns values sampled from the probability distribution of the mean selection coefficient for a given set of bed regions.\n" +
			"Selection calculated from sequence variability from an input VCF.\n" +
			"Uses derived allele frequencies from VCFs with ancestor allele annotations by default.\n" +
			"Implements the Metropolis-Hastings algorithm to evaluate a Markov Chain Monte Carlo approximation of the selection coefficient distribution.\n" +
			"Usage:\n" +
			"selectionMCMC input.vcf out.txt\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var iterations *int = flag.Int("iterations", 100, "Number of MCMC iterations.")
	var muZero *float64 = flag.Float64("muZero", 0, "Starting position for the mean selection coefficient parameter mu.")
	var sigmaZero *float64 = flag.Float64("sigmaZero", 0.1, "Starting value for the selection coefficient distribution variance parameter sigma.")
	var muStep *float64 = flag.Float64("muStep", 0.2, "Step size for the mean selection coefficient parameter mu.")
	var sigmaStep *float64 = flag.Float64("sigmaStep", 0.1, "Step size for the selection coefficient variance parameter sigma.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var unPolarized *bool = flag.Bool("unPolarized", false, "Disable the requirement for ancestor annotation and use unpolarized site frequency spectrum. Use with caution.")
	var divergenceAscertainment *bool = flag.Bool("divergenceAscertainment", false, "Make a divergence-based ascertainment correction.")
	var fixedSigma *bool = flag.Bool("fixedSigma", false, "When true, the selection coefficient variance parameter sigma stays fixed at sigmaZero and is not treated as an independent hyperparameter.")
	var integralError *float64 = flag.Float64("integralError", 1e-7, "Set the error threshold for numerical integration.")
	var verbose *int = flag.Int("verbose", 0, "Set to 1 or 2 to reveal different levels of debug print statements to standard output.")
	var sigmaPriorAlpha *float64 = flag.Float64("sigmaPriorAlpha", 2, "Sets the alpha parameter for the Gamma-distributed prior distribution for the hyperparameter sigma.")
	var sigmaPriorBeta *float64 = flag.Float64("sigmaPriorBeta", 10, "Sets the beta parameter for the Gamma-distributed prior distribution for the hyperparameter sigma.")
	var muPriorMean *float64 = flag.Float64("muPriorMean", 0, "Sets the mean of the normally-distributed prior distribution for the hyperparameter mu.")
	var muPriorSigma *float64 = flag.Float64("muPriorSigma", 3, "Sets the standard deviation of he normally-distributed prior distribution for the hyperparameter mu.")
	var includeRef *bool = flag.Bool("includeRef", false, "Include the reference genome allele state as a datapoint in the allele frequency spectrum.")

	flag.Usage = usage
	//log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	log.SetFlags(0)
	flag.Parse()

	s := popgen.McmcSettings{
		Iterations:              *iterations,
		MuStep:                  *muStep,
		MuZero:                  *muZero,
		SigmaStep:               *sigmaStep,
		SigmaZero:               *sigmaZero,
		SetSeed:                 *setSeed,
		UnPolarized:             *unPolarized,
		DivergenceAscertainment: *divergenceAscertainment,
		FixedSigma:              *fixedSigma,
		D:                       1, //D is hardcoded as 1 for now. This represents the size of the ascertainment subset.
		IntegralError:           *integralError,
		Verbose:                 *verbose,
		SigmaPriorAlpha:         *sigmaPriorAlpha,
		SigmaPriorBeta:          *sigmaPriorBeta,
		MuPriorMean:             *muPriorMean,
		MuPriorSigma:            *muPriorSigma,
		IncludeRef:              *includeRef,
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	vcfFile := flag.Arg(0)
	outFile := flag.Arg(1)
	selectionMcmc(vcfFile, outFile, s)
}
