// Command Group: "Sequence Evolution & Reconstruction"

// Generates a newick tree file from an input dot format tree
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/tree"
	"log"
)

func dotToNewick(inFile string, outFile string, verbose bool) {
	t := tree.ParseDot(inFile, verbose)
	err := tree.WriteNewick(outFile, t)
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"dotToNewick - Generates a newick tree file from an input dot format tree.\n" +
			"Usage:\n" +
			"dotToNewick input.dot output.nh\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var verbose *bool = flag.Bool("verbose", false, "Enables debug prints.")
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outFile := flag.Arg(1)

	dotToNewick(infile, outFile, *verbose)
}
