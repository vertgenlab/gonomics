// Command Group: "Sorting"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/genomeGraph"
)

func usage() {
	fmt.Print(
		"sortGraph - Topologically sorts nodes in a genome graph (.gg) file. \n" +
			"Usage:\n" +
			" sortGraph input.gg output.gg\n")
	flag.PrintDefaults()
}

func sortGraph(inFile string, outFile string) {
	graph := genomeGraph.Read(inFile)
	graph = genomeGraph.SortGraph(graph)
	genomeGraph.Write(outFile, graph)
}

func main() {
	expectedNumArgs := 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("ERROR: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	sortGraph(inFile, outFile)
}
