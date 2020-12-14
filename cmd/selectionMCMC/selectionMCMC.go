package main

import (
	"flag"
	"fmt"
	"log"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/common"
)

func selectionMCMC(filename string, outFile string, muZero float64, sigmaZero float64, iterations int, randSeed bool, setSeed int64) {
	common.RngSeed(randSeed, setSeed)
	data := popgen.VcfToAFS(filename)
	popgen.MetropolisHastings(data, muZero, sigmaZero, iterations, outFile)
}

func usage() {
	fmt.Print(
		"selectionMCMC - Returns values sampled from the probability distribution of the mean selection coefficient for a given set of bed regions.\n" +
			"Selection calculated from sequence variability from an input gVCF.\n" +
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
	var randSeed *bool = flag.Bool("randSeed", false, "Uses a random seed for the RNG.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")

	flag.Usage = usage
	//log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	log.SetFlags(0)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	vcfFile := flag.Arg(0)
	outFile := flag.Arg(1)
	selectionMCMC(vcfFile, outFile, *muZero, *sigmaZero, *iterations, *randSeed, *setSeed)
}
