package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func faNoGap(inFile string, outFile string) {
	records := fasta.Read(inFile)

	beds := fasta.UngappedRegionsAll(records)

	bed.Write(outFile, beds, 3)
}

func usage() {
	fmt.Print(
		"faNoGap - bed file containing regions outside gaps\n" +
			"Usage:\n" +
			" faNoGap input.fa output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	faNoGap(inFile, outFile)
}
