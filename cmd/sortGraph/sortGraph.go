package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
)

func usage() {
	fmt.Print(
		"sortGraph - Topologically sorts nodes in a genome graph (.gg) file. \n" +
			"Usage:\n" +
			" sortGraph input.gg output.gg\n")
	flag.PrintDefaults()
}

func sortGraph(inFile string, outFile string) {
	graph := simpleGraph.Read(inFile)
	graph = simpleGraph.SortGraph(graph)
	simpleGraph.Write(outFile, graph)
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
