// Command Group: "Data Simulation"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/genomeGraph"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
)

func girafSimulate(graph *genomeGraph.GenomeGraph, numReads int, readLen int, randSeed int64, numSomaticSNV int, AlleleFrequency float64, outFile string, outputSam bool) {
	reads := genomeGraph.RandGiraf(graph, numReads, readLen, randSeed)
	var samReads []sam.Sam
	if numSomaticSNV != 0 {
		genomeGraph.RandSomaticMutations(graph, reads, numSomaticSNV, AlleleFrequency, randSeed)
	}

	if outputSam == true {
		for i := 0; i < len(reads); i++ {
			samReads = append(samReads, genomeGraph.GirafToSam(reads[i]))
		}
		//samHeader := genomeGraph.NodesHeader(graph.Nodes)
		sam.Write(outFile, samReads, sam.Header{})
	} else {
		giraf.Write(outFile, reads)
	}
}

func usage() {
	fmt.Print(
		"girafSimulate - Returns a file of giraf alignments for a input genome graph.\n" +
			"Usage:\n" +
			" girafSimulate [options] input.gg \n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var numReads *int = flag.Int("numReads", 100, "Number of giraf alignments to simulate.")
	var readLen *int = flag.Int("readLen", 150, "Length of each alignment.")
	var seed *int64 = flag.Int64("seed", 0, "Seed used to generate alignments, default is time-based seed")
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
	graph := genomeGraph.Read(inFile)

	//TODO add remove block once GirafToSam is complete
	if *outputSam == true {
		log.Fatalln("ERROR: Sam output is still in development")
	}

	girafSimulate(graph, *numReads, *readLen, *seed, *numSomaticSNV, *AlleleFrequency, outFile, *outputSam)
}
