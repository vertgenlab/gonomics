// Command Group: "Data Conversion"

// Reference-based diploid assembly of aligned short reads
package main

import (
	"flag"
	"fmt"
	"log"
)

const bufferSize = 10_000_000

func usage() {
	fmt.Print(
		"samAssembler - Reference-based diploid assembly of aligned short reads.\n" +
			"Can be used in three modes:\n\t'build' generates diploid assemblies.\n" +
			"\t'score' validates assembly accuracy with a five way alignment including the known divergent sequences\n" +
			"\t'prior' constructs an empirical prior for output diploid genotypes based on a maximum likelihood estimate from the input reads.\n" +
			"Enter: 'samAssembler build' OR 'samAssembler score' OR 'samAssembler prior' to view usage and options.\n")
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Expecting at least one argument, received 0.")
	}

	var mode = flag.Arg(0)

	switch mode {
	case "build":
		parseBuildArgs()
	case "score":
		parseScoreArgs()
	case "prior":
		parsePriorArgs()
	default:
		log.Fatalf("Unknown mode. samAssembler can be run with the first argument as 'build', 'prior', or 'score'.")
	}
}
