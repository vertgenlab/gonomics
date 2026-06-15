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
    var 

	for i := 0; i < len(b); i++ {
		fasta.(alnFile, outFile, b[i].ChromStart, b[i].ChromEnd, noMask, lineLength, false) //hard code endOfAlignment as false, as we are getting end positions from the beds.
	}
}

func usage() {
	fmt.Print(
		"multFaBedDivergence - Calculates the divergence (number of mutations) between two sequences in an alignment.\n" +
			"All bed entries must be on the same chromosome to interface with multiFa file.\n" +
			"The number of divergences will be added to the Score field of the bed file.\n" +
			"Usage:\n" +
			"multiFaVisualizeBeds in.bed aln.mfa out.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var noMask *bool = flag.Bool("noMask", false, "Converts all bases to upper case.")
	var outFormat *bool = flag.Bool("outFormatName", false, "Uses the name column as the outfile name (name.txt).")
	var lineLength *int = flag.Int("lineLength", 100, "Sets to length of each alignment line.")
	var outDir *string = flag.String("outDir", "", "Set a path for the output files. Should end with \"/\".")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	bedFile := flag.Arg(0)
	alnFile := flag.Arg(1)

	multFaVisualizeBeds(bedFile, alnFile, *outFormat, *noMask, *lineLength, *outDir)
}
