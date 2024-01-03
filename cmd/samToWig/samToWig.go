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

type Settings struct {
	SamFileName        string
	ChromSizesFileName string
	OutFileName        string
	FragLength         int
	DefaultValue       float64
	Deletions          bool
}

func samToWig(s Settings) {
	ref := chromInfo.ReadToMap(s.ChromSizesFileName)
	inChan, _ := sam.GoReadToChan(s.SamFileName)
	answer := wig.MakeSkeleton(ref, s.DefaultValue)
	var j int
	var currentBed bed.Bed
	var currBeds []bed.Bed

	for currSam := range inChan {
		if s.FragLength != -1 {
			currentBed = convert.SamToBedFrag(currSam, s.FragLength, ref)
			if currentBed.Chrom != "" {
				convert.BedReadUpdateWig(answer, currentBed)
			}
		} else if s.Deletions {
			currBeds = convert.SamToBedWithDeletions(currSam)
			for j = range currBeds {
				if currBeds[j].Chrom != "" {
					convert.BedReadUpdateWig(answer, currBeds[j])
				}
			}
		} else {
			currentBed = convert.SamToBed(currSam)
			if currentBed.Chrom != "" {
				convert.BedReadUpdateWig(answer, currentBed)
			}
		}
	}
	wig.Write(s.OutFileName, answer)
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
	var defaultValue *float64 = flag.Float64("defaultValue", 0, "Specify a default value for the output wig in positions without read coverage.")
	var deletions *bool = flag.Bool("deletions", false, "Won't fill in wig values over deletions. Sam records that have deletions will contribute 0 to the deletion locations."+
		"Not compatible with fragLength.")

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

	s := Settings{
		SamFileName:        inFile,
		ChromSizesFileName: reference,
		OutFileName:        outFile,
		FragLength:         *fragLength,
		DefaultValue:       *defaultValue,
		Deletions:          *deletions,
	}

	samToWig(s)
}
