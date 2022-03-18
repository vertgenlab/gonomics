// Command Group: "Data Simulation"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	//"github.com/vertgenlab/gonomics/exception"
	//"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
	"math/rand"
)

func simulateAFS(popSize uint, mutRate int, genTime uint, genomeSize int, outFile string, setSeed int64) {
	rand.Seed(setSeed)
	c := simulate.GoSimulateAFS(popSize, mutRate, genTime, genomeSize)
	//var err error

	fasta.Write(outFile, c)

	

	fmt.Println(c)
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
	var popSize *uint = flag.Uint("N", 10, "Specifies the population size in the simulation.")
	var mutRate *int = flag.Int("m", 1, "Specifies the genome-wide mutation rate in the simulation.")
	var genTime *uint = flag.Uint("t", 10, "Specifies the generation time in the simulation")
	var genomeSize *int = flag.Int("g", 100, "Specifies the genome size in base-pair in the simulation")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	outFile := flag.Arg(0)

	simulateAFS(*popSize, *mutRate, *genTime, *genomeSize, outFile, *setSeed)
}
