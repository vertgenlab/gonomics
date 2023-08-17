// Command Group: "Data Conversion"

// Convert the output of the aidenlab's straw command to a bedpe format
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/hic"
	"log"
)

func strawToBedpe(strawFile string, outFile string, chrom string, binSize int, interChrom string) {
	var thisRec bedpe.BedPe
	var err error

	straw := hic.GoReadToChan(strawFile)
	out := fileio.EasyCreate(outFile)

	if interChrom == "" {
		for s := range straw {
			thisRec = bedpe.BedPe{A: bed.Bed{Chrom: chrom, ChromStart: s.Bin1Start, ChromEnd: s.Bin1Start + binSize, Score: s.ContactScore, FieldsInitialized: 8}, B: bed.Bed{Chrom: chrom, ChromStart: s.Bin2Start, ChromEnd: s.Bin2Start + binSize, Score: s.ContactScore, FieldsInitialized: 8}}
			log.Printf("this record: %v", thisRec.B.Chrom)
			bedpe.WriteToFileHandle(out, thisRec)
		}
	} else {
		for s := range straw {
			thisRec = bedpe.BedPe{A: bed.Bed{Chrom: chrom, ChromStart: s.Bin1Start, ChromEnd: s.Bin1Start + binSize, Score: s.ContactScore, FieldsInitialized: 8}, B: bed.Bed{Chrom: interChrom, ChromStart: s.Bin2Start, ChromEnd: s.Bin2Start + binSize, Score: s.ContactScore, FieldsInitialized: 8}}
			bedpe.WriteToFileHandle(out, thisRec)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"strawToBedpe - convert the output of the aidenlab's straw command to a bedpe format." +
			"\n" +
			"Usage:\n" +
			"strawToBedpe [options] file.straw out.bedpe chrom\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var binSize *int = flag.Int("binSize", 5000, "A binSize must be provided for the resolution of the straw file output if it is not 5000.")
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
