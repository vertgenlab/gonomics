package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/sam"
	//"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func totalAlignedBases(filename string) int {
	samFile := fileio.EasyOpen(filename)
	defer samFile.Close()
	var done bool = false
	var aln sam.Aln
	var alignedBases int

	sam.ReadHeader(samFile)

	for aln, done = sam.ReadNext(samFile); done != true; aln, done = sam.ReadNext(samFile) {
		if aln.Cigar[0].Op != '*' {
			alignedBases += cigar.MatchLength(aln.Cigar)
		}
	}
	return alignedBases
}

func samCoverage(samFileName string, noGapFileName string, outFile string) {
	noGap := bed.Read(noGapFileName)
	genomeSize := bed.TotalSize(noGap)
	alignedBases := totalAlignedBases(samFileName)

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
