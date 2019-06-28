package main 

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/convert"
	"log"
)

func samToBed(infile string, outfile string, paired *bool, fragLength *int64) {
	fmt.Printf("Paired: %b\n", paired)
	fmt.Printf("fragLength: %v\n", fragLength)

	if *paired && *fragLength != int64(-1) {
		log.Fatalf("Invalid entry. Cannot be both paired and have a fixed frag size.")
	}

	records, err := sam.Read(infile)
	if err != nil {
		log.Fatal(err)
	}
	var outBed []*bed.Bed

/* TODO: Write paired command 
	if paired {
		outBed = convert.SamToBedPaired(records)
	} else*/
	if *fragLength != int64(-1) {
		outBed = convert.SamToBedFrag(records, *fragLength)
	} else {
		outBed = convert.SamToBed(records)
	}

	bed.Write(outfile, outBed, 4)
}

func usage() {
	fmt.Print(
		"samToBed - Converts sam to bed\n" +
		"Usage:\n" +
		" samToBed input.sam output.bed\n" +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
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
	outFile := flag.Arg(1)

	samToBed(inFile, outFile, paired, fragLength)
}