// Command Group: "Sequence Evolution & Reconstruction"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/reconstruct"
	"log"
)

func chimpAncestorRecon(inFile string, outFile string, messyToN bool) {
	records := fasta.Read(inFile)
	output := append(records, reconstruct.ChimpAncestorRecon(records, messyToN))
	fasta.Write(outFile, output)
}

func usage() {
	fmt.Print(
		"chimpAncestorRecon - Like quickPrimateRecon, but returns a chimp-biased HCA estimate from a four-way primate multiFa alignment.\n" +
			"Must be in the following order: Human, Bonobo, Chimp, Orangutan, Gorilla\n" +
			"Usage:\n" +
			"  quickPrimateRecon input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var messyToN *bool = flag.Bool("messyToN", false, "SEts messy bases to Ns in the output file.")
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

	chimpAncestorRecon(inFile, outFile, *messyToN)
}
