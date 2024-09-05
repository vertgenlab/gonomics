// Command Group: "FASTA and Multi-FASTA Tools"

// Removes all columns in a multi fasta alignment that are not variable
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func mfaReduce(inFilename, outFilename, bedFilename, chrom string) {
	aln := fasta.Read(inFilename)
	var answer []fasta.Fasta
	var answerBedPos []int      // answerBedPos is []int not [][]bed.Bed to avoid circular dependencies, aka trying to import bed package in fasta package
	var answerBedNames []string // answerBedName holds the referenceSpecies1base_alignSpecies2base of each segregating site
	if bedFilename != "" {
		var answerBed []bed.Bed
		var currentBed bed.Bed
		var chromStartAlnPos, chromStartRefPos int
		lastAlnPosConverted := 0
		lastRefPosConverted := 0
		answer, answerBedPos, answerBedNames = fasta.SegregatingSitesWithBed(aln)
		for i := 0; i < len(answerBedPos); i++ {
			chromStartAlnPos = answerBedPos[i]
			chromStartRefPos = fasta.AlnPosToRefPosCounter(aln[0], chromStartAlnPos, lastRefPosConverted, lastAlnPosConverted)
			lastAlnPosConverted = chromStartAlnPos
			lastRefPosConverted = chromStartRefPos
			currentBed = bed.Bed{Chrom: chrom, ChromStart: chromStartRefPos, ChromEnd: chromStartRefPos + 1, Name: answerBedNames[i], Score: chromStartAlnPos, FieldsInitialized: 5} // Name field is reference species name, Score field is AlnPos
			answerBed = append(answerBed, currentBed)
		}
		bed.Write(bedFilename, answerBed)
	} else {
		answer = fasta.SegregatingSites(aln)
	}
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
	var bedFilename *string = flag.String("bedFilename", "", "Output the positions of variable sites in the reference species into a bed file with this name. Variable sites are reported 1 base/line")
	var chrom *string = flag.String("chrom", "", "Required when using -bedFilename, to specify the chromosome name of the reference species in the output bed file")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	if *bedFilename != "" && *chrom == "" {
		flag.Usage()
		log.Fatalf("Error: using -bedFilename without -chrom\n")
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	mfaReduce(inFile, outFile, *bedFilename, *chrom)
}
