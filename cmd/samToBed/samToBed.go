package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

func samToBed(infile string, reference string, outfile string, paired bool, fragLength int64) {
	log.Printf("Paired: %b\n", paired)
	log.Printf("fragLength: %v\n", fragLength)

	if paired && fragLength != int64(-1) {
		log.Fatalf("Invalid entry. Cannot be both paired and have a fixed frag size.")
	}
	ref := chromInfo.ReadToMap(reference)
	records, err := sam.Read(infile)
	common.ExitIfError(err)
	var outBed []*bed.Bed

	/* TODO: Write paired command
	if paired {
		outBed = convert.SamToBedPaired(records)
	} else*/
	if fragLength != int64(-1) {
		outBed = convert.SamToBedFrag(records, fragLength, ref)
	} else {
		outBed = convert.SamToBed(records)
	}

	bed.Write(outfile, outBed, 4)
}

func usage() {
	fmt.Print(
		"samToBed - Converts sam to bed\n" +
			"Usage:\n" +
			" samToBed input.sam reference.chrom.sizes output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var paired *bool = flag.Bool("windowSize", false, "Specifies paired end reads")
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

	samToBed(inFile, reference, outFile, *paired, *fragLength)
}
