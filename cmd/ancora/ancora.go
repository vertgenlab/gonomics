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
		"ANCoRA - Ancient-DNA Nucleotide-damage Correction and Reference-guided Assembly.\n" +
			"This program builds reference-guided assemblies from input short-read sequencing libraries.\n" +
			"Optimized for ancient DNA libraries, and learns ancient DNA damage profiles for genotype correction.\n" +
			"Can be used in three modes:\n\t'build' generates diploid assemblies.\n" +
			"\t'score' validates assembly accuracy with a five way alignment including the known divergent sequences\n" +
			"\t'prior' constructs an empirical prior for output diploid genotypes based on a maximum likelihood estimate from the input reads.\n" +
			"Enter: 'ancora build' OR 'ancora score' OR 'ancora prior' to view usage and options.\n")
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
		log.Fatalf("Unknown mode. ANCoRA can be run with the first argument as 'build', 'prior', or 'score'.")
	}
}
