// Command Group: "FASTA and Multi-FASTA Tools"

// Provides the mutational distance between two sequences human-readable multiple alignments for all entries in a bed file
package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/bed"
)

func multFaBedDivergence(bedFile string, alnFile string, outFile string) {
	b := bed.Read(bedFile)

	for i := 0; i < len(b); i++ {
		fasta.(alnFile, outFile, b[i].ChromStart, b[i].ChromEnd, noMask, lineLength, false)
	}
}

func usage() {
	fmt.Print(
		"multFaBedDivergence - Calculates the divergence (number of mutations) between two sequences in an alignment.\n" +
			"All bed entries must be on the same chromosome to interface with the multiFa file.\n" +
			"The number of divergences will appear in the Score field of the bed file.\n" +
			"Usage:\n" +
			"multiFaBedDivergence in.bed aln.mfa out.bed\n" +
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
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	bedFilename := flag.Arg(0)
	alnFilename := flag.Arg(1)
	outBedFilename := flag.Arg(2)

	multiFaBedDivergence(bedFilename, alnFilename, outBedFilename)
}
