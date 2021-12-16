// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/exception"
	"log"
)

func samToWig(samFileName string, reference string, outfile string, fragLength int) {
	log.Printf("fragLength: %d\n", fragLength)

	ref := chromInfo.ReadToMap(reference)

	samFile := fileio.EasyOpen(samFileName)
	var err error

	var done bool = false
	sam.ReadHeader(samFile)
	var outBed []bed.Bed
	var aln sam.Sam

	var outWig []wig.Wig
	var currentBed bed.Bed

	if fragLength != -1 {
		for aln, done = sam.ReadNext(samFile); done != true; aln, done = sam.ReadNext(samFile) {
			currentBed = convert.SamToBedFrag(aln, fragLength, ref)
			if currentBed.Chrom != "" {
				outBed = append(outBed, currentBed)
			}
		}
	} else {
		for aln, done = sam.ReadNext(samFile); done != true; aln, done = sam.ReadNext(samFile) {
			currentBed = convert.SamToBed(aln)
			if currentBed.Chrom != "" {
				outBed = append(outBed, currentBed)
			}
		}
	}
	err = samFile.Close()
	exception.PanicOnErr(err)

	outWig = convert.BedReadsToWig(outBed, ref)
	log.Printf("Length of outWig: %d", len(outWig))
	wig.Write(outfile, outWig)
}

func usage() {
	fmt.Print(
		"samToWig - Converts sam to wig\n" +
			"Usage:\n" +
			" samToWig input.sam reference.chrom.sizes output.wig\n" +
			" Currently fills in Wig values over deletions.\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var fragLength *int = flag.Int("fragLength", -1, "Specifies the fragment length for ChIP-Seq, must be greater than or equal to read length")

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

	samToWig(inFile, reference, outFile, *fragLength)
}
