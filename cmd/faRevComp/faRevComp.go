package main

import (
	"flag"
	"fmt"
	"github.com/craiglowe/gonomics/fasta"
	"log"
)

func faRevComp(inFile string, outFile string, keepName bool) {
	records, err := fasta.Read(inFile)
	if err != nil {
		log.Fatal(err)
	}

	fasta.ReverseComplementAll(records)

	if keepName == false {
		fasta.AppendToNameAll(records, "_revComp")
	}

	fasta.Write(outFile, records)
}

func usage() {
	fmt.Print(
		"faRevComp - reverse complement the sequences in a fasta file\n" +
			"Usage:\n" +
			" faRevComp input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var keepName *bool = flag.Bool("keepName", false, "keep the sequence name the same instead of appending \"_revComp\"")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	faRevComp(inFile, outFile, *keepName)
}
