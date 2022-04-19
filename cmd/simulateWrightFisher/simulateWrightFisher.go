// Command Group: "Data Simulation"

package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/simulate"
)

func simulateWrightFisher(outFile string, set popgen.WrightFisherSettings) {
	rand.Seed(set.SetSeed)
	wf := simulate.SimulateWrightFisher(set)
	if set.Fasta {
		popgen.WriteTSV(outFile, wf)
	} else {
		fasta.Write(outFile, wf.Fasta)
	}

}

func usage() {
	fmt.Print(
		"simulateWrightFisher - Returns a file of fasta from a simulation of given population size, mutation rate, and generation time.\n" +
			"Usage:\n" +
			" simulateWrightFisher output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	var popSize *int = flag.Int("N", 100, "Specifies the population size in the simulation.")
	var mutRate *float64 = flag.Float64("m", 1e-1, "Specifies the genome-wide mutation rate in the simulation.")
	var numGen *int = flag.Int("t", 100, "Specifies the number of generations passed in the simulation")
	var genomeSize *int = flag.Int("g", 1, "Specifies the genome size in base-pair in the simulation")
	var rFitness *float64 = flag.Float64("w", 2, "Specifies the relative fitness of the derived allele over ancestral allele. Must be greater or equal than zero")
	var gcContent *float64 = flag.Float64("gc", 0.5, "Specifies the GC content for the simulated ancestral sequence")

	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var verbose *bool = flag.Bool("verbose", false, "Verbose: Use this flag for debugging purposes")
	var fasta *bool = flag.Bool("f", false, "Use this flag to set output to be a multialignment fasta")
	var vcf *bool = flag.Bool("v", false, "Use this flag to set output to be a vcf")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	s := popgen.WrightFisherSettings{
		PopSize:    *popSize,
		MutRate:    *mutRate,
		NumGen:     *numGen,
		GenomeSize: *genomeSize,
		RFitness:   *rFitness,
		GcContent:  *gcContent,
		SetSeed:    *setSeed,
		Verbose:    *verbose,
		Fasta:      *fasta,
		Vcf:        *vcf,
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	outFile := flag.Arg(0)

	simulateWrightFisher(outFile, s)
}
