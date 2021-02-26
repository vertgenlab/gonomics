package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
)

func simulateBed(regionCount int, simLength int, noGapFile string, outFile string, randSeed bool, setSeed int64) {
	common.RngSeed(randSeed, setSeed)
	noGap := bed.Read(noGapFile)
	c := simulate.GoSimulateBed(noGap, regionCount, simLength)
	out := fileio.EasyCreate(outFile)
	defer out.Close()

	for i := range c {
		bed.WriteBed(out.File, &i, 5)
	}
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
	var randSeed *bool = flag.Bool("randSeed", false, "Uses a random seed for the RNG.")
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

	simulateBed(*regionCount, *Length, inFile, outFile, *randSeed, *setSeed)
}
