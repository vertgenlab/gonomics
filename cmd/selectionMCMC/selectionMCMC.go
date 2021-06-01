// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/popgen"
	"log"
)

func selectionMCMC(filename string, outFile string, s popgen.McmcSettings) {
	common.RngSeed(s.RandSeed, s.SetSeed)
	data, err := popgen.VcfToAfs(filename, !s.UnPolarized) //VcfToAFS is written with polarized as the argument for clarity, so the bool is flipped here.
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
	var sigmaZero *float64 = flag.Float64("sigmaZero", 1, "Starting value for the selection coefficient distribution variance parameter sigma.")
	var muStep *float64 = flag.Float64("muStep", 0.2, "Step size for the mean selection coefficient parameter mu.")
	var sigmaStep *float64 = flag.Float64("sigmaStep", 50, "Step size for the selection coefficient variance parameter sigma.")
	var randSeed *bool = flag.Bool("randSeed", false, "Uses a random seed for the RNG.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var unPolarized *bool = flag.Bool("unPolarized", false, "Disable the requirement for ancestor annotation and use unpolarized site frequency spectrum. Use with caution.")
	var derived *bool = flag.Bool("derived", false, "Make a divergence-based ascertainment correction for regions enriched for derived alleles (i.e. HAQERs, HARs, or other fast-evolving regions).")
	var ancestral *bool = flag.Bool("ancestral", false, "Make a divergence-based ascertainment correction for regions enriched for ancestral alleles (i.e. UCEs or other highly conserved regions).")
	var integralError *float64 = flag.Float64("integralError", 1e-7, "Set the error threshold for numerical integration.")

	flag.Usage = usage
	//log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	log.SetFlags(0)
	flag.Parse()

	options := popgen.McmcSettings{
		Iterations:    *iterations,
		MuStep:        *muStep,
		MuZero:        *muZero,
		SigmaStep:     *sigmaStep,
		SigmaZero:     *sigmaZero,
		RandSeed:      *randSeed,
		SetSeed:       *setSeed,
		UnPolarized:   *unPolarized,
		Derived:       *derived,
		Ancestral:     *ancestral,
		IntegralError: *integralError,
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	vcfFile := flag.Arg(0)
	outFile := flag.Arg(1)
	selectionMCMC(vcfFile, outFile, options)
}
