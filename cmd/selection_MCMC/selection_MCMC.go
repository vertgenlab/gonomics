package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"log"
)

func selectionMCMC(filename string, outfile string, muZero float64, sigmaZero float64, iterations int, logSpace bool) {
	data := popgen.GVCFToAFS(filename)
	muList, sigmaList, acceptList := popgen.MetropolisHastings(data, muZero, sigmaZero, iterations, logSpace)
	out := fileio.EasyCreate(outfile)
	defer out.Close()
	fmt.Fprintf(out, "Iteration\tMu\tSigma\tAccept\n")
	for i := 0; i < len(muList); i++ {
		fmt.Fprintf(out, "%v\t%e\t%e\t%t\n", i, muList[i], sigmaList[i], acceptList[i])
	}
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
	var logSpace *bool = flag.Bool("logSpace", false, "This argument converts internal calculations to logSpace. Recommended for large datasets.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	vcfFile := flag.Arg(0)
	outFile := flag.Arg(1)
	selectionMCMC(vcfFile, outFile, *muZero, *sigmaZero, *iterations, *logSpace)
}
