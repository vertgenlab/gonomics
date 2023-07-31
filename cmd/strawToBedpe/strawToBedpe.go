// Command Group: "Data Conversion"

// Convert the output of the aidenlab's straw command to a bedpe format
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/hic"
	"log"
)

func strawToBedpe(strawFile string, outFile string, chrom string, binSize int, interChrom string) {
	straw := hic.Read(strawFile)
	var out = make([]bedpe.BedPe, len(straw))

	if interChrom == "" {
		for s := range straw {
			out[s] = bedpe.BedPe{bed.Bed{Chrom: chrom, ChromStart: straw[s].Bin1Start, ChromEnd: straw[s].Bin1Start + binSize, Score: straw[s].ContactScore}, bed.Bed{Chrom: chrom, ChromStart: straw[s].Bin2Start, ChromEnd: straw[s].Bin2Start + binSize, Score: straw[s].ContactScore}}
		}
	} else {
		for s := range straw {
			out[s] = bedpe.BedPe{bed.Bed{Chrom: chrom, ChromStart: straw[s].Bin1Start, ChromEnd: straw[s].Bin1Start + binSize, Score: straw[s].ContactScore}, bed.Bed{Chrom: interChrom, ChromStart: straw[s].Bin2Start, ChromEnd: straw[s].Bin2Start + binSize, Score: straw[s].ContactScore}}
		}
	}

	bedpe.Write(outFile, out)
}

func usage() {
	fmt.Print(
		"strawToBedpe - convert the output of the aidenlab's straw command to a bedpe format." +
			"\n" +
			"Usage:\n" +
			"\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var binSize *int = flag.Int("binsize", 5000, "A binSize must be provided for the resolution of the straw file output if it is not 5000.")
	var interChrom *string = flag.String("interChrom", "", "If the straw file contains regions from two different chromosomes, this option should hold the chromosome for the second bin's chromosome name. The argument ''chrom'' should hold the first bin's chromosome name.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	strawFile := flag.Arg(0)
	outFile := flag.Arg(1)
	chrom := flag.Arg(3)

	strawToBedpe(strawFile, outFile, chrom, *binSize, *interChrom)
}
