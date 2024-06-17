// Command Group: "WIG Tools"

// a collection of tools for manipulating Wig files.
package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Printf(
		"wigTools - a collection of tools for manipulating Wig files.\n" +
			"Usage:\n" +
			"\twigTools peaks in.wig out.bed\n" +
			"\twigTools filter in.wig genome.chrom.sizes out.wig\n" +
			"\tOR\n" +
			"\twigTools toTrainingSet in.wig genome.fa train.txt validate.txt test.txt\n" +
			"\tOR\n" +
			"\twigTools math in.wig chrom.sizes out.file\n" +
			"\tOR\n" +
			"\twigTools stats in.wig chrom.sizes noGap.bed output.tsv\n" +
			"Enter a subcommand to view options.\n")
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Error: user must specify a wigTools subcommand.\n")
	}

	switch flag.Arg(0) {
	case "filter":
		parseFilterArgs()
	case "peaks":
		parsePeaksArgs()
	case "toTrainingSet":
		parseToTrainingSetArgs()
	case "math":
		parseMathArgs()
	case "stats":
		parseStatsArgs()
	default:
		flag.Usage()
		log.Fatalf("Error: unrecognized subcommand: %v.\n", flag.Arg(0))
	}
}
