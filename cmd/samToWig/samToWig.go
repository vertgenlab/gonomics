// Command Group: "Data Conversion"

// Converts sam or bam to wig
package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/wig"
)

func samToWig(samFileName string, reference string, outfile string, fragLength int, deletions bool) {
	log.Printf("fragLength: %d\n", fragLength)

	ref := chromInfo.ReadToMap(reference)

	inChan, _ := sam.GoReadToChan(samFileName)

	var outBed, currBeds []bed.Bed

	var outWig []wig.Wig
	var currentBed bed.Bed

	for i := range inChan {
		if fragLength != -1 {
			currentBed = convert.SamToBedFrag(i, fragLength, ref)
			if currentBed.Chrom != "" {
				outBed = append(outBed, currentBed)
			}
		} else if deletions {
			currBeds = convert.SamToBedWithDeletions(i)
			for j := range currBeds {
				if currBeds[j].Chrom != "" {
					outBed = append(outBed, currBeds[j])
				}
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
	var deletions *bool = flag.Bool("deletions", false, "Won't fill in wig values over deletions. Sam records that have deletions will contribute 0 to the deletion locations."+
		"Not compatible with fragLength")

	flag.Usage = usage

	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	if *fragLength != -1 && *deletions {
		log.Fatalln("ERROR: -fragLength is not compatible with -deletions")
	}

	inFile := flag.Arg(0)
	reference := flag.Arg(1)
	outFile := flag.Arg(2)

	samToWig(inFile, reference, outFile, *fragLength, *deletions)
}
