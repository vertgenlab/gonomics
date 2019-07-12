package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/wig"
	"log"
)

func samToWig(infile string, reference string, outfile string, paired bool, fragLength int64) {
	fmt.Printf("Paired: %t\n", paired)
	fmt.Printf("fragLength: %d\n", fragLength)
	if paired && fragLength != int64(-1) {
		log.Fatalf("Invalid entry. Cannot be both paired and have a fixed frag size.")
	}

	records, err := sam.Read(infile)
	if err != nil {
		log.Fatal(err)
	}

	ref := chromInfo.ReadToMap(reference)

	var outBed []*bed.Bed
	var outWig []*wig.Wig

	if fragLength != int64(-1) {
		outBed = convert.SamToBedFrag(records, fragLength, ref)
	} else {
		outBed = convert.SamToBed(records)
	}

	outWig = convert.BedToWig(outBed, ref)
	fmt.Printf("Length of outWig: %d", len(outWig))
	wig.Write(outfile, outWig)
}

func usage() {
	fmt.Print(
		"samToWig - Converts sam to wig\n" +
			"Usage:\n" +
			" samToWig input.sam reference.chrom.sizes output.wig\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var paired *bool = flag.Bool("paired", false, "Specifies paired end reads")
	var fragLength *int64 = flag.Int64("fragLength", -1, "Specifies the fragment length for ChIP-Seq")

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

	samToWig(inFile, reference, outFile, *paired, *fragLength)
}
