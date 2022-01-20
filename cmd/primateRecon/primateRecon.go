// Command Group: "Sequence Evolution & Reconstruction"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/reconstruct"
	"log"
)

func primateRecon(infile string, outfile string, messyToN bool) {
	records := fasta.Read(infile)
	output := append(records, reconstruct.PrimateRecon(records, messyToN))
	fasta.Write(outfile, output)
}

func usage() {
	fmt.Print(
		"primateRecon - Returns maximum likelihood sequence from an HBCGO primate alignment\n" +
			"Usage:\n" +
			"primateRecon input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var messyToN *bool = flag.Bool("messyToN", false, "Sets messy bases to Ns in the output file.")
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

	primateRecon(inFile, outFile, *messyToN)
}
