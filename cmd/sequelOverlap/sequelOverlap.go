package main

import (
	"flag"
	"fmt"
	"log"
)

func usage() {
	fmt.Print(
		"sequelOverlap - A tool to find non/overlapping genomic regions\n\n" +
			"Usage:\n" +
			"  sequelOverlap [options] in.file out.file\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

func main() {

	//flags:
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
	} else {
		//inFile, outFile := flag.Arg(0), flag.Arg(1)
		//overlapSelect()
	}
}

/*
Considering how many options we might have I think we should make this separate from the main function
func overlapSelect(flag string) {
	switch flag {
	case:

	case:
	}
}*/
