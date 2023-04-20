// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/fasta"
)

func assemblyStats(infile string, outfile string, countLowerAsGaps bool) {
	N50, halfGenome, genomeLength, largestContig, numContigs := fasta.AssemblyStats(infile, countLowerAsGaps)
	fasta.WriteAssemblyStats(infile, outfile, N50, halfGenome, genomeLength, largestContig, numContigs)
}

func usage() {
	fmt.Print(
		"assemblyStats - Provides information about the number of scaffolds, including the N50, number of scaffolds, and distribution of lengths of assembled scaffolds.\n" +
			"Usage:\n" +
			"assemblyStats input.fa output.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var countLower *bool = flag.Bool("countLowerAsGaps", false, "Lower case letters count as gaps and break contigs.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	infile := flag.Arg(0)
	outfile := flag.Arg(1)
	assemblyStats(infile, outfile, *countLower)
}
