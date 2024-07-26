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

func simulateVcf(s SimulateVcfSettings) {
	rand.Seed(s.SetSeed)
	simulate.VcfToFile(s.Alpha, s.NumAlleles, s.NumSites, s.OutFile, s.BoundAlpha, s.BoundBeta, s.BoundMultiplier, s.RefFile, s.HasRef)
}

func usage() {
	fmt.Print(
		"simulateVcf - Contains functions for simulating VCF data.\n" +
			"Usage:\n" +
			" simulateBed [options] output.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

type SimulateVcfSettings struct {
	OutFile         string
	Alpha           float64
	NumAlleles      int
	NumSites        int
	SetSeed         int64
	BoundAlpha      float64
	BoundBeta       float64
	BoundMultiplier float64
	RefFile         *string
	HasRef          bool
}

func main() {
	var expectedNumArgs int = 1
	SimulateVcfFlags := flag.NewFlagSet("SimulateVcf", flag.ExitOnError)
	var numSites *int = SimulateVcfFlags.Int("numSites", 10, "Specifies the number of simulated variants.")
	var setSeed *int64 = SimulateVcfFlags.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var alpha *float64 = SimulateVcfFlags.Float64("alpha", 0.01, "Specifies the selection parameter alpha for drawing individual gVCF alleles from a stationarity distribution.")
	var numAlleles *int = SimulateVcfFlags.Int("numAlleles", 0, "Specifies the number of alleles for gVCF samples.")
	var boundAlpha *float64 = SimulateVcfFlags.Float64("boundAlpha", 0.001, "Set the alpha parameter for the bounding function.")
	var boundBeta *float64 = SimulateVcfFlags.Float64("boundBeta", 0.001, "Set the beta parameter for the bounding function.")
	var boundMultiplier *float64 = SimulateVcfFlags.Float64("boundMultiplier", 10000, "Set the multiplier for the bounding function.")
	var refFile *string = SimulateVcfFlags.String("refFile", "", "Specify a reference Fasta file.")
	var hasRef bool = false
	SimulateVcfFlags.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	SimulateVcfFlags.Parse()

	if strings.ToLower(refFile) != "" {
		hasRef = true
	}
	if len(SimulateVcfFlags.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(SimulateVcfFlags.Args()))
	}
	outFile := SimulateVcfFlags.Arg(0)

	s := SimulateVcfSettings{
		OutFile:         outFile,
		Alpha:           *alpha,
		NumAlleles:      *numAlleles,
		NumSites:        *numSites,
		SetSeed:         *setSeed,
		BoundAlpha:      *boundAlpha,
		BoundBeta:       *boundBeta,
		BoundMultiplier: *boundMultiplier,
		RefFile:         *refFile,
		HasRef:          *hasRef,
	}

	simulateVcf(s)
}
