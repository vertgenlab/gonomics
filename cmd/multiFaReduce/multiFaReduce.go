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

func mfaReduce(inFilename, outFilename, bedFilename, chrom string, refStart int) {
	aln := fasta.Read(inFilename)
	var answer []fasta.Fasta
	if bedFilename != "" {
		var answerBed []bed.Bed
		answer, answerBed = bed.SegregatingSites(aln, chrom, refStart)
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
	var bedFilename *string = flag.String("bedFilename", "", "Output the positions of variable sites into a bed file with this name. Positions are reported on the reference species (the top/first sequence in the multiFa alignment). Variable sites are reported 1 base/line.")
	var chrom *string = flag.String("chrom", "", "Required when using -bedFilename, to specify the chromosome name of the reference species in the output bed file.")
	var refStart *int = flag.Int("refStart", 0, "Optional when using -befFilename, to set the reference position for the beginning of the input multiFa alignment.")

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

	mfaReduce(inFile, outFile, *bedFilename, *chrom, *refStart)
}
