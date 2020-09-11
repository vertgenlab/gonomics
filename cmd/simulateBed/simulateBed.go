package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"log"
	"math/rand"
	"time"
)

func simulateBed(regionCount int, simLength int64, noGapFile string, outFile string, randSeed bool, setSeed int64) {
	if randSeed && setSeed != -1 {
		log.Fatalf("Cannot use a set seed and also a random seed.")
	}
	if randSeed {
		rand.Seed(time.Now().UnixNano())
	} else if setSeed != -1 {
		rand.Seed(setSeed)
	}

	noGap := bed.Read(noGapFile)
	answer := make([]*bed.Bed, 0)
	var Length, tmp, chromWindows int64
	var totalWindows int

	//count total viable windows
	for i := 0; i < len(noGap); i++ {
		Length = noGap[i].ChromEnd - noGap[i].ChromStart

		if Length >= simLength {
			totalWindows = totalWindows + int(Length-simLength)
		}
	}
	//fmt.Printf("totalWindows: %d\n", totalWindows)

	for i := 0; i < regionCount; i++ {
		tmp = int64(common.RandIntInRange(0, totalWindows))
		//fmt.Printf("Random number is %d\n", tmp)
		for j := 0; j < len(noGap); j++ {
			Length = noGap[j].ChromEnd - noGap[j].ChromStart
			chromWindows = Length - simLength + 1
			//is chrom big enough?
			if chromWindows < 1 {
				break
			}
			if tmp-chromWindows > 0 {
				tmp = tmp - chromWindows
			} else {
				//DEBUG: fmt.Printf("Got one\n")
				answer = append(answer, &bed.Bed{Chrom: noGap[j].Chrom, ChromStart: noGap[j].ChromStart + tmp - 1, ChromEnd: noGap[j].ChromStart + tmp - 1 + simLength, Name: noGap[j].Name})
				break
			}
		}
	}
	bed.Write(outFile, answer, 4)
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
	var Length *int64 = flag.Int64("L", 1000, "Specifies the length of simulated regions.")
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
