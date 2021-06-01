// Command Group: "Data Simulation"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
)

func simulateVcf(alpha float64, n int, k int, outFile string, randSeed bool, setSeed int64) {
	common.RngSeed(randSeed, setSeed)
	simulate.SimulateVcf(alpha, n, k, outFile)
}

func usage() {
	fmt.Print(
		"simulateVcf - Contains functions for simulating VCF data.\n" +
			"Usage:\n" +
			" simulateBed [options] output.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	var numSites *int = flag.Int("numSites", 10, "Specifies the number of simulated bed regions.")
	var randSeed *bool = flag.Bool("randSeed", false, "Uses a random seed for the RNG.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var alpha *float64 = flag.Float64("alpha", 0.01, "Specifies the selection parameter alpha for drawing individual gVCF alleles from a stationarity distribution.")
	var numAlleles *int = flag.Int("numAlleles", 0, "Specifies the number of alleles for gVCF samples.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	outFile := flag.Arg(0)

	simulateVcf(*alpha, *numAlleles, *numSites, outFile, *randSeed, *setSeed)
}
