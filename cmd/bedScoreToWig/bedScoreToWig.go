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

func bedScoreToWig(infile string, reference string, outfile string) {
	rec := bed.Read(infile)
	ref := chromInfo.ReadToMap(reference)

	outWig := convert.BedScoreToWig(rec, ref)
	fmt.Println("Writing the outfile.")
	wig.Write(outfile, outWig)
}

func usage() {
	fmt.Print(
		"bedScoreToWig - Converts bed score to wig\n" +
			"Usage:\n" +
			"bedScoreToWig input.bed reference.chrom.sizes output.wig\n" +
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

	bedScoreToWig(inFile, reference, outFile)
}