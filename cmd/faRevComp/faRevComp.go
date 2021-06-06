// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func faRevComp(inFile string, outFile string, keepName bool) {
	records := fasta.Read(inFile)

	fasta.ReverseComplementAll(records)

	if keepName == false {
		for i := range records {
			records[i].Name = records[i].Name + "_revComp"
		}
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
