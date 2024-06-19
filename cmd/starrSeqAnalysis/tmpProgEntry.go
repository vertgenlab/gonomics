package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Print("cellrangerBam -- Takes in a cellranger bam file of STARR-seq reads and parses the extra flags field to pull out the " +
		"representative read for each UMI and which construct it belongs to. Multiple GEM wells from the same STARR-seq experiment can be provided in a comma-delimited list " +
		"in the 'inFile' field. The output is a tab-delimited table of read-counts for each constructs.\n" +
		"NOTE: The default behavior of this function works best with STARR-seq libraries where constructs don't have much similarity with each other.\n" +
		"For libraries that need barcoding (like GWAS or cross-species comparisons) use the -bed option with a bed file corresponding to barcode regions.\n" +
		"Usage: \n" +
		"cellrangerBam [options] inFile outFile\n\n")
	flag.PrintDefaults()
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Error: user must specify a starrSeqAnalysis subcommand.\n")
	}

	switch flag.Arg(0) {
	case "makeRef":
		parseMakeRefArgs()
	case "inputSeq":
		parseInputSeqArgs()
	case "outputSeq":
		parseOutputSeqArgs()
	case "makeOligoPool":
		parseMakeLibraryArgs()
	default:
		flag.Usage()
		log.Fatalf("Error: unrecognized subcommand: %v\n", flag.Arg(0))
	}
}
