package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/bed"
	//"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func samCoverage(samFileName string, noGapFileName string, outFile string) {
	noGap := bed.Read(noGapFileName)
	genomeSize := bed.TotalSize(noGap)
	alignedBases := sam.TotalAlignedBases(samFileName)

	coverage := (float64(alignedBases) / float64(genomeSize))
	out := fileio.EasyCreate(outFile)
	defer out.Close()

	fmt.Fprintf(out, "Aligned Bases: %d\n", alignedBases)
	fmt.Fprintf(out, "Genome Length: %d\n", genomeSize)
	fmt.Fprintf(out, "Coverage: %v\n", coverage)
}

func usage() {
	fmt.Print(
		"samCoverage - Calculates genome coverage as the quotient of aligned bases in a sequencing dataset to the total length of ungapped genomic regions in the reference genome.\n" +
		"Usage:\n" +
		"samCoverage input.sam nogap.bed outfile.txt\n" +
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

	samFileName := flag.Arg(0)
	noGapFileName := flag.Arg(1)
	outFile := flag.Arg(2)


	samCoverage(samFileName, noGapFileName, outFile)
}
