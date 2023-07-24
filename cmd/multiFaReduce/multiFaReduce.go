// Command Group: "FASTA and Multi-FASTA Tools"

// Removes all columns in a multi fasta alignment that are not variable
package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/fasta"
)

func mfaReduce(inFilename, outFilename string) {
	aln := fasta.Read(inFilename)
	answer := fasta.SegregatingSites(aln)
	fasta.Write(outFilename, answer)
}

func usage() {
	fmt.Print(
		"mfaReduce - mfaReduce removes all columns in a multi fasta alignment that are not variable\n" +
			"Usage:\n" +
			" mfaReduce input.mfa output.mfa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	mfaReduce(inFile, outFile)
}
