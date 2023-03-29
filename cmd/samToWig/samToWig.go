// Command Group: "Data Conversion"

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

func samToWig(samFileName string, reference string, outfile string, fragLength int) {
	log.Printf("fragLength: %d\n", fragLength)

	ref := chromInfo.ReadToMap(reference)

	inChan, _ := sam.GoReadToChan(samFileName)

	var outBed []bed.Bed

	var outWig []wig.Wig
	var currentBed bed.Bed

	for i := range inChan {
		if fragLength != -1 {
			currentBed = convert.SamToBedFrag(i, fragLength, ref)
			if currentBed.Chrom != "" {
				outBed = append(outBed, currentBed)
			}
		} else {
			currentBed = convert.SamToBed(i)
			if currentBed.Chrom != "" {
				outBed = append(outBed, currentBed)
			}
		}
	}

	outWig = convert.BedReadsToWig(outBed, ref)
	//log.Printf("Length of outWig: %d", len(outWig))
	wig.Write(outfile, outWig)
}

func usage() {
	fmt.Print(
		"samToWig - Converts sam or bam to wig\n" +
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
