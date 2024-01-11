// Command Group: "Sequence Evolution & Reconstruction"

// Simulate sequence evolution on an input fasta from the root of a newick tree to the leaf nodes
package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Print(
		"simulateEvol - a suite of tools for molecular evolution simulation.\n" +
			"Select a subcommand to view options.\n" +
			"Usage:\n" +
			"\tsimulateEvol genic tree.newick in.fasta out.fasta\n" +
			"\tOR\n" +
			"\tsimulateEvol withIndels seq.fasta outFile.fasta\n" +
			"\tOR\n" +
			"\tsimulateEvol nonCoding out.fasta\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Error: user must specify a simulateEvol subcommand.\n")
	}

	switch flag.Arg(0) {
	case "genic":
		parseGenicArgs()
	case "withIndels":
		parseWithIndelsFlags()
	case "nonCoding":
		parseNonCodingArgs()
	default:
		flag.Usage()
		log.Fatalf("Unrecognized subcommand: %v.\n", flag.Arg(0))
	}
}
