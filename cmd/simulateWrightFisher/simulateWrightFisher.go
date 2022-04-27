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
		fasta.Write(outFile, wf.Fasta)
	} else {
		popgen.WriteTSV(outFile, wf)
	}

}

func usage() {
	fmt.Print(
		"simulateWrightFisher - simulate a multiallelic, haplotic Wright-Fisher population (discrete, non-overlapping generations).\n" +
			"\t\tBy default, returns a tsv file containing allele frequencies of all alleles at all sites in all generation.\n" +
			"\t\tUse -f or -v flag to switch the output to fasta or vcf, respectively\n" +
			"Usage:\n" +
			"\tsimulateWrightFisher [OPTIONS] output.tsv\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	var popSize *int = flag.Int("N", 100, "Specifies the population size in the simulation. (Default 100)")
	var mutRate *float64 = flag.Float64("m", 1e-1, "Specifies the genome-wide, uniform mutation rate in the simulation. (Default 0.1)")
	var numGen *int = flag.Int("t", 500, "Specifies the number of generations passed in the simulation. (Defailt 500)")
	var genomeSize *int = flag.Int("g", 1, "Specifies the genome size in base-pair in the simulation (Default 1)")
	var rFitness *float64 = flag.Float64("w", 1, "Specifies the relative fitness of the derived allele over ancestral allele. Must be greater or equal than zero. (Default 1)")
	var gcContent *float64 = flag.Float64("gc", 0.5, "Specifies the GC content for the simulated ancestral sequence. (Default 0.5)")
	var initFreq *string = flag.String("i", "", "Specifies the initial frequencies for all alleles for all sites, as well as ancestral allele.\n"+
		"Accepts comma-separated string of frequencies of A,C,G,T respectively, and then an ancestral allele.\n"+
		"Example: -i 0.25,0.25,0.25,0.25,A gives frequencies of A,C,G,T to be 0.25,0.25,0.25,0.25, and set A as the ancestral allele for all sites.\n"+
		"This overrides -gc.")
	var fitnessString *string = flag.String("W", "", "Specifies relative frequencies of each allele (A,C,G,T), and specifies ancestral allele.\n"+
		"Same format as -i. Example: -W 1,1.02,0.98,1,A. The ancestral allele MUST BE THE SAME as input from -i, if applicable. This overrides -w.")

	var setSeed *int64 = flag.Int64("setSeed", 1, "Use a specific seed for the RNG. (Default 1)")
	var verbose *bool = flag.Bool("verbose", false, "Verbose: Use this flag for debugging purposes")
	var fasta *bool = flag.Bool("f", false, "Use this flag to set output to be a multialignment fasta")
	var vcf *bool = flag.Bool("v", false, "Use this flag to set output to be a vcf")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	s := popgen.WrightFisherSettings{
		PopSize:       *popSize,
		MutRate:       *mutRate,
		NumGen:        *numGen,
		GenomeSize:    *genomeSize,
		RFitness:      *rFitness,
		GcContent:     *gcContent,
		InitFreq:      *initFreq,
		FitnessString: *fitnessString,
		SetSeed:       *setSeed,
		Verbose:       *verbose,
		Fasta:         *fasta,
		Vcf:           *vcf,
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	outFile := flag.Arg(0)

	simulateWrightFisher(outFile, s)
}
