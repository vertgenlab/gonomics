package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Printf(
		"samInfo - a collection of tools for extracting information from SAM/BAM files." +
			"Usage:\n" +
			"\tsamInfo readLength in.sam/bam out.txt\n" +
			"\tOR\n" +
			"\tsamInfo coverage in.sam/bam histogram.txt statSummary.txt\n" +
			"\tOR\n" +
			"\tsamInfo mapq in.sam/bam histogram.txt\n" +
			"Enter a subcommand to view options.\n")
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Error: user must specify a samInfo subcommand.\n")
	}

	switch flag.Arg(0) {
	case "readLength":
		parseReadLengthArgs()
	case "coverage":
		parseCoverageArgs()
	case "mapq":
		parseMapqArgs()
	default:
		flag.Usage()
		log.Fatalf("Error: unrecognized subcommand: %v.\n", flag.Arg(0))
	}
}
