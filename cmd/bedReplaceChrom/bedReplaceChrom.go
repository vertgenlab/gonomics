// Command Group: "BED Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"log"
)

func bedReplaceChrom(infile string, outfile string, chromName string) {
	var records []*bed.Bed = bed.Read(infile)

	for i := 0; i < len(records); i++ {
		records[i].Chrom = chromName
	}
	bed.Write(outfile, records, 5)
}

func usage() {
	fmt.Print(
		"bedReplaceChrom - replaces the Name field of a bed file with input.\n" +
			"Usage:\n" +
			"bedFilter input.bed output.bed chromName\n" +
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

	infile := flag.Arg(0)
	outfile := flag.Arg(1)
	chromName := flag.Arg(2)

	bedReplaceChrom(infile, outfile, chromName)
}
