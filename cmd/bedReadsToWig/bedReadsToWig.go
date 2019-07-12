package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/wig"
	"log"
)

func bedReadsToWig(infile string, reference string, outfile string) {
	rec := bed.Read(infile)
	ref := chromInfo.ReadToMap(reference)

	outWig := convert.BedReadsToWig(rec, ref)
	wig.Write(outfile, outWig)
}

func usage() {
	fmt.Print(
		"bedReadsToWig - Converts bed reads to wig\n" +
			"Usage:\n" +
			"bedReadsToWig input.bed reference.chrom.sizes output.wig\n" +
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
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	reference := flag.Arg(1)
	outFile := flag.Arg(2)

	bedReadsToWig(inFile, reference, outFile)
}
