package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/sort"
	"log"
)

func usage() {
	fmt.Print(
		"sortGiraf - External sort of giraf records based on topological ordering of nodes in input graph.\n" +
			"Usage:\n" +
			" sortGiraf -o output.giraf -g graph.gg [options] input.giraf \n" +
			"options:\n")
	flag.PrintDefaults()
}

func sortGiraf(girafFile string, graphFile string, linesPerChunk int, outFile string) {
	graph := simpleGraph.Read(graphFile)
	topoOrder := simpleGraph.GetSortOrder(graph)
	sort.ExternalMergeSort(girafFile, topoOrder, linesPerChunk, outFile)
}

func main() {
	var expectedNumArgs = 1
	var outFile *string = flag.String("o", "", "Write output to a file [.giraf].")
	var graph *string = flag.String("g", "", "Graph reference used for generating input giraf.")
	var linesPerChunk *int = flag.Int("chunkSize", 1000, "Number of giraf records to use for each tmp file written to disk.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("ERROR: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	if *graph == "" || *outFile == "" {
		flag.Usage()
		log.Fatalf("ERROR: Must include parameters for both -o and -g")
	}

	inFile := flag.Arg(0)

	sortGiraf(inFile, *graph, *linesPerChunk, *outFile)
}
