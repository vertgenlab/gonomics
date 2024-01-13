package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Printf(
		"oboTools - a collection of tools for manipulating files in the Open Biomedical Ontologies (OBO) file format.\n" +
			"Subcommands include:\n" +
			"\tmapping - a tool for creating tabular maps between Obo fields.\n" +
			"Enter subcommands to view specific usage statements:\n" +
			"\toboTools mapping\n")
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Error: User must specify an oboTools subcommand to view usage statements.\n")
	}

	switch flag.Arg(0) {
	case "mapping":
		parseMappingArgs()
	default:
		log.Fatalf("Error: unrecognized subcommand: %v.\n", flag.Arg(0))
	}
}
