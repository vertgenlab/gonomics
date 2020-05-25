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
			" sortGiraf -o output.giraf -g graph.gg [options] input.giraf graphReference.gg output.giraf \n" +
			"options:\n")
	flag.PrintDefaults()
}

func sortGiraf(girafFile string, graphFile string, linesPerChunk int, outFile string) {
	graph := simpleGraph.Read(graphFile)
	topoOrder := simpleGraph.GetSortOrder(graph)
	sort.ExternalMergeSort(girafFile, topoOrder, linesPerChunk, outFile)
}

func main() {
	var expectedNumArgs = 3
	//TODO: add flag that allows user to define custom prefix for tmp files
	var linesPerChunk *int = flag.Int("chunkSize", 1000000, "Number of giraf records to use for each tmp file written to disk.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("ERROR: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	graph := flag.Arg(1)
	outFile := flag.Arg(2)

	sortGiraf(inFile, graph, *linesPerChunk, outFile)
}
