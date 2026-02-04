// Command Group: "Data Simulation"

// Contains functions for simulating VCF data
package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"
	"strings"

	"github.com/vertgenlab/gonomics/simulate"
)

func simulateVcf(s Settings) {
	rand.Seed(s.SetSeed)
	simulate.VcfToFile(s.Alpha, s.NumAlleles, s.NumSites, s.OutFile, s.BoundAlpha, s.BoundBeta, s.BoundMultiplier, s.RefFile, s.HasRef)
}

func usage() {
	fmt.Print(
		"simulateVcf - Simulate VCF format variant data.\n" +
			"Sample genotype allele frequencies are drawn from a derived allele frequency spectrum\n" +
			"specified by a selection parameter alpha (alpha = 2N_{e}s).\n" +
			"By default, all variants are (A -> T) at sequential positions.\n" +
			"Alternatively, the user can simulate variants from a real reference genome using the 'RefFile' option.\n" +
			"Usage:\n" +
			" simulateVcf [options] output.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	OutFile         string
	Alpha           float64
	NumAlleles      int
	NumSites        int
	SetSeed         int64
	BoundAlpha      float64
	BoundBeta       float64
	BoundMultiplier float64
	RefFile         string
	HasRef          bool
}

func main() {
	var expectedNumArgs int = 1
	var numSites *int = flag.Int("numSites", 10, "Specifies the number of simulated variants.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var alpha *float64 = flag.Float64("alpha", 0.01, "Specifies the selection parameter alpha for drawing individual gVCF alleles from a stationarity distribution.")
	var numAlleles *int = flag.Int("numAlleles", 10, "Specifies the number of alleles for gVCF samples.")
	var boundAlpha *float64 = flag.Float64("boundAlpha", 0.001, "Set the alpha parameter for the bounding function.")
	var boundBeta *float64 = flag.Float64("boundBeta", 0.001, "Set the beta parameter for the bounding function.")
	var boundMultiplier *float64 = flag.Float64("boundMultiplier", 10000, "Set the multiplier for the bounding function.")
	var refFile *string = flag.String("refFile", "", "Specify a reference Fasta file.")
	var hasRef bool = false

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *refFile != "" {
		hasRef = true
	}
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	outFile := flag.Arg(0)

	s := Settings{
		OutFile:         outFile,
		Alpha:           *alpha,
		NumAlleles:      *numAlleles,
		NumSites:        *numSites,
		SetSeed:         *setSeed,
		BoundAlpha:      *boundAlpha,
		BoundBeta:       *boundBeta,
		BoundMultiplier: *boundMultiplier,
		RefFile:         *refFile,
		HasRef:          hasRef,
	}

	simulateVcf(s)
}
