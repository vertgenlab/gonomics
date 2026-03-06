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
			"\tpFaTools extractBed in.pfa regions.bed out.pFa\n" +
			"\tOR\n" +
			"\tpFaTools sample in.pfa outDir\n" +
			"\tOR\n" +
			"\tPFaTools visualize in.pfa start end out.txt\n" +
			"\tOR\n" +
			"\tPFaTools faToPfa inFile out.pfa InputType\n" +
			"\tOR\n" +
			"\tPFaTools vcfToPfa inFile.vcf refFile.fa out.pfa\n" +
			"\tOR\n" +
			"\tPFaTools distTrack infile1.pfa infile2.pfa/fa out.wig\n" +
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
	case "faToPfa":
		parseFaToPfaArgs()
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
	case "vcfToPfa":
		parseVcfToPfaArgs()
	case "distTrack":
		parseDistTrackArgs()
	default:
		flag.Usage()
		log.Fatalf("Error: unrecognized subcommand: %v.\n", flag.Arg(0))
	}
}
