// Command Group: "Sequence Evolution & Reconstruction"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/tree"
)

func dotToNewick(inFile string, outFile string, verbose bool) {
	t := tree.ParseDot(inFile, verbose)
	secondErr := tree.WriteNewick(outFile, t)
	common.ExitIfError(secondErr)
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
