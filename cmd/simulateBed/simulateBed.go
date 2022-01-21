// Command Group: "Data Simulation"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
	"math/rand"
)

func simulateBed(regionCount int, simLength int, noGapFile string, outFile string, setSeed int64) {
	rand.Seed(setSeed)
	noGap := bed.Read(noGapFile)
	c := simulate.GoSimulateBed(noGap, regionCount, simLength)
	out := fileio.EasyCreate(outFile)
	var err error

	for i := range c {
		bed.WriteBed(out.File, i)
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"simulateBed - Returns a file of random bed regions of an input bed file.\n" +
			"Usage:\n" +
			" simulateBed input.bed output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var Length *int = flag.Int("L", 1000, "Specifies the length of simulated regions.")
	var regionCount *int = flag.Int("N", 10, "Specifies the number of simulated bed regions.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	simulateBed(*regionCount, *Length, inFile, outFile, *setSeed)
}
