// Command Group: "Data Simulation"

package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simulate"
)

func simulateWrightFisher(popSize int, mutRate float64, numGen int, genomeSize int, rFitness float64, gcContent float64, outFile string, setSeed int64, verbose int) {
	rand.Seed(setSeed)
	c := simulate.SimulateWrightFisher(popSize, mutRate, numGen, genomeSize, rFitness, gcContent, verbose)
	fasta.Write(outFile, c)
}

func usage() {
	fmt.Print(
		"simulateAFS - Returns a file of fasta from a simulation of given population size, mutation rate, and generation time.\n" +
			"Usage:\n" +
			" simulateAFS output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	var popSize *int = flag.Int("N", 5, "Specifies the population size in the simulation.")
	var mutRate *float64 = flag.Float64("m", 5e-1, "Specifies the genome-wide mutation rate in the simulation.")
	var numGen *int = flag.Int("t", 10, "Specifies the number of generations passed in the simulation")
	var genomeSize *int = flag.Int("g", 10, "Specifies the genome size in base-pair in the simulation")
	var rFitness *float64 = flag.Float64("w", 1, "Specifies the relative fitness of the derived allele over ancestral allele")
	var gcContent *float64 = flag.Float64("gc", 0.5, "Specifies the GC content for the simulated ancestral sequence")

	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var verbose *int = flag.Int("verbose", 0, "Set to 1 to reveal statements to stdout for debugging")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	outFile := flag.Arg(0)

	simulateWrightFisher(*popSize, *mutRate, *numGen, *genomeSize, *rFitness, *gcContent, outFile, *setSeed, *verbose)
}
