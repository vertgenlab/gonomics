package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
	"time"
)

func simulateGiraf(graph *simpleGraph.SimpleGraph, numReads int, readLen int, randSeed int64, numSomaticSNV int, AlleleFrequency float64, outFile string, outputSam bool) {
	reads := simpleGraph.RandGiraf(graph, numReads, readLen, randSeed)
	var samReads []*sam.SamAln

	if numSomaticSNV != 0 {
		simpleGraph.RandSomaticMutations(graph, reads, numSomaticSNV, AlleleFrequency, randSeed)
	}

	if outputSam == true {
		for i := 0 ; i < len(reads) ; i++ {
			samReads = append(samReads, simpleGraph.GirafToSam(reads[i]))
		}
		samHeader := simpleGraph.NodesHeader(graph.Nodes)
		sam.Write(outFile, &sam.Sam{samHeader, samReads})
	} else {
		giraf.Write(outFile, reads)
	}
}

func usage() {
	fmt.Print(
		"simulateGiraf - Returns a file of giraf alignments for a input genome graph.\n" +
			"Usage:\n" +
			" simulateGiraf [options] input.gg \n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var time = time.Now().UnixNano()
	var expectedNumArgs int = 2
	var numReads *int = flag.Int("numReads", 100, "Number of giraf alignments to simulate.")
	var readLen *int = flag.Int("readLen", 150, "Length of each alignment.")
	var seed *int64 = flag.Int64("seed", time, "Seed used to generate alignments, default is time-based seed")
	var numSomaticSNV *int = flag.Int("somaticSNV", 0, "Number of simulated somatic SNVs.")
	var AlleleFrequency *float64 = flag.Float64("somaticAF", 0.2, "Approx. allele frequency of simulated somatic mutations.")
	var outputSam *bool = flag.Bool("outputSam", false, "Output file in sam format")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("ERROR: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)
	graph := simpleGraph.Read(inFile)

	simulateGiraf(graph, *numReads, *readLen, *seed, *numSomaticSNV, *AlleleFrequency, outFile, *outputSam)

}