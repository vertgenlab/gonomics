package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/ontology/obo"
	"log"
)

func oboToDot(oboFile string, term string, output string) {
	oboRecords, _ := obo.Read(oboFile, false)

	obo.SubtreeToDot(output, term, oboRecords)
}

func usage() {
	fmt.Print(
		"oboToDot takes in an obo file and the id belonging to a GO term in the obo file and returns the dot " +
			"formatted subtree belonging to that ID " +
			"Usage:\n" +
			"oboToDot input.obo GoIdName output.dot\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	oboFile := flag.Arg(0)
	id := flag.Arg(1)
	out := flag.Arg(2)

	oboToDot(oboFile, id, out)
}
