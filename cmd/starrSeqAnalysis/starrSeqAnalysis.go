package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Print("starrSeqAnalysis -- A suite of tools to analysis single cell STARR-seq data\n\n" +
		"Subcommands:\n" +
		"makeRef\n" +
		"inputSeq\n" +
		"outputSeq\n\n" +
		"type a subcommand to print the help message for that command:\n" +
		"starrSeqAnyalsis [subcommand]\n\n")
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
