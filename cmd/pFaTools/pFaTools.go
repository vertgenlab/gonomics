// Command Group: "FASTA and Multi-FASTA Tools"

// A collection of tools for manipulating pFasta files
package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Printf(
		"pFaTools - a collection of tools for manipulating pFasta files.\n" +
			"Usage:\n" +
			"\tpFaTools entropyTrack in.pfa out.wig\n" +
			"\tOR\n" +
			"\tpFaTools extract in.pfa chrom start end out.pFa\n" +
			"\tOR\n" +
			"\tpFaTools extractBed in.pfa regions.bed outDir.pFa\n" +
			"\tOR\n" +
			"\tpFaTools sample in.pfa outDir.Fa\n" +
			"\tOR\n" +
			"\tPFaTools visualize in.pfa start end outDir.txt\n" +
			"Enter a subcommand to view options.\n")
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Error: user must specify a pFaTools subcommand.\n")
	}

	switch flag.Arg(0) {
	case "entropyTrack":
		parseEntropyTrackArgs()
	case "extract":
		parseExtractArgs()
	case "extractBed":
		parseExtractBedArgs()
	case "sample":
		parseSampleArgs()
	case "visualize":
		parseVisualizeArgs()
	default:
		flag.Usage()
		log.Fatalf("Error: unrecognized subcommand: %v.\n", flag.Arg(0))
	}
}
