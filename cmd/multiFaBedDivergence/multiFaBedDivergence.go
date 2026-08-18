// Command Group: "FASTA and Multi-FASTA Tools"

// Provides the mutational distance between two sequences human-readable multiple alignments for all entries in a bed file
package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/bed"
        "github.com/vertgenlab/gonomics/fasta"
)

func multiFaBedDivergence(bedFile string, alnFile string, outFile string) {
	beds := bed.Read(bedFile)
	aln := fasta.Read(alnFile)

	var refIdx, alnIdx int = 0, 0
	var chromStartAln, chromEndAln, divergence int

	for i := 0; i < len(beds); i++ {
		chromStartAln = fasta.RefPosToAlnPosCounter(aln[0], beds[i].ChromStart, refIdx, alnIdx)
		refIdx = beds[i].ChromStart
		alnIdx = chromStartAln
		chromEndAln = fasta.RefPosToAlnPosCounter(aln[0], beds[i].ChromEnd, refIdx, alnIdx)
		refIdx = beds[i].ChromEnd
		alnIdx = chromEndAln
		divergence = fasta.PairwiseMutationDistanceInRange(aln[1], aln[2], chromStartAln, chromEndAln)
		beds[i].Score = divergence
		if beds[i].FieldsInitialized < 4 {
			beds[i].Name = fmt.Sprintf("element%d", i)
			beds[i].FieldsInitialized = 4
		}
		if beds[i].FieldsInitialized < 5 {
			beds[i].FieldsInitialized = 5
		}
		log.Printf("%s\n", beds[i])
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
	var expectedNumArgs int = 3
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
