package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Printf(
		"gtfTools - a collection of tools for Gene Transfer Format files.\n" +
			"Usage:\n" +
			"\tgtfTools filter in.gtf out.gtf\n" +
			"\tgtfTools toBed in.gtf out.bed\n" +
			"Enter a subcommand to view options.\n")
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Error: user must specify a gtfTools subcommand.\n")
	}

	switch flag.Arg(0) {
	case "filter":
		parseFilterArgs()
	case "toBed":
		parseToBedArgs()
	default:
		flag.Usage()
		log.Fatalf("Error: unrecognized subcommand: %v.\n", flag.Arg(0))
	}
}
