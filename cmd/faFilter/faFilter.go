package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func faFilter(infile string, outfile string, name string) {
	records := fasta.Read(infile)
	var outlist []*fasta.Fasta

	for i := 0; i < len(records); i++ {
		if records[i].Name == name {
			outlist = append(outlist, records[i])
		}
	}
	fasta.Write(outfile, outlist)
}

func usage() {
	fmt.Print(
		"faFilter - Returns a filtered fasta with headers matching a given name.\n" +
			"Usage:\n" +
			" faFindFast input.fa output.fa name\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)
	name := flag.Arg(2)

	faFilter(inFile, outFile, name)
}
