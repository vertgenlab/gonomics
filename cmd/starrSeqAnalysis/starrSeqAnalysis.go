package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Print("starrSeqAnalysis -- A suite of tools to analysis single cell STARR-seq data\n\n" +
		"\tstarrSeqAnalysis makeRef up.fa constructs.fa down.fa outfilePrefix\n" +
		"\tstarrSeqAnalysis inputSeq [options] in.sam in.bed out.txt\n" +
		"\tstarrSeqAnalysis outputSeq [options] in.sam out.txt\n" +
		"enter a subcommand to view options and help message:\n")
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
