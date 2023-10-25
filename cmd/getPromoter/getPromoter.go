// Command Group: "GTF Tools"

// getPromoter will take a list of unique genes of interest and return a bed file of the promoter region proceeding
// the position of the TSS for each isoform by the amount specified with upstream and following the TSS by the amount specified by downstream.
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/fileio"
	gtf2 "github.com/vertgenlab/gonomics/gtf"
	"log"
)

func getPromoter(genes string, info string, outBed string, chrom string, upstream int, downstream int) {
	geneNames := fileio.Read(genes)
	geneInfo := gtf2.Read(info)
	chr := chromInfo.ReadToMap(chrom)

	answer := gtf2.FindPromoter(geneNames, upstream, downstream, geneInfo, chr)

	bed.Write(outBed, answer)
}

func usage() {
	fmt.Print(
		"getPromoter will take a list of unique genes of interest and return a bed file of the promoter region proceeding\n" +
			"the position of the TSS for each isoform in a given gtf file by the amount specified with upstream and the amount following the TSS as specified by downstream. " +
			"This program is strand aware and requires a chrom sizes file for the genome of interest.\n" +
			"Usage:\n" +
			"getPromoter [options] uniqueGenes.txt geneInfo.gtf out.bed chrom.sizes\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4

	var upstream *int = flag.Int("upstream", 1000, "Specifies the amount of bases to return in the bed entry before the TSS")
	var downstream *int = flag.Int("downstream", 200, "Specifies the amount of bases to return in the bed entry after the TSS")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	uniqueGenes := flag.Arg(0)
	geneInfo := flag.Arg(1)
	outFile := flag.Arg(2)
	chromInfoName := flag.Arg(3)

	getPromoter(uniqueGenes, geneInfo, outFile, chromInfoName, *upstream, *downstream)
}
