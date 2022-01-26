// Command Group: "Sorting"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/genomeGraph"
	"github.com/vertgenlab/gonomics/sort"
	"log"
)

func usage() {
	fmt.Print(
		"girafSort - External sort of giraf records based on topological ordering of nodes in input graph.\n" +
			"Usage:\n" +
			" girafSort -o output.giraf -g graph.gg [options] input.giraf graphReference.gg output.giraf \n" +
			"options:\n")
	flag.PrintDefaults()
}

func girafSort(girafFile string, graphFile string, linesPerChunk int, outFile string) []uint32 {
	graph := genomeGraph.Read(graphFile)
	topoOrder := genomeGraph.GetSortOrder(graph)
	sort.GirafExternalMergeSort(girafFile, topoOrder, linesPerChunk, outFile)
	return topoOrder // topo order is not deterministic so we have a return for testing purposes only
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

	girafSort(inFile, graph, *linesPerChunk, outFile)
}
