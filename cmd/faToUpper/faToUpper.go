// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func faToUpper(inFile string, outFile string, chrom string) {
	records := fasta.Read(inFile)

	if chrom != "" {
		for i := 0; i < len(records); i++ {
			if records[i].Name == chrom {
				fasta.ToUpper(records[i])
			}
		}
	} else {
		for i := 0; i < len(records); i++ {
			fasta.ToUpper(records[i])
		}
	}
	fasta.Write(outFile, records)
}

func usage() {
	fmt.Print(
		"faToUpper - reformat a fasta file in upper case\n" +
			"Usage:\n" +
			" faToUpper input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var chrom *string = flag.String("chrom", "", "Only moves entries matching a specific name to uppercase.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	faToUpper(inFile, outFile, *chrom)
}
