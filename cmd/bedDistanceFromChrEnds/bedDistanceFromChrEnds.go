// Command Group: "BED Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

func bedDistanceFromChrEnds(inFile string, chromFile string, outFile string) {
	records := bed.Read(inFile)
	ref := chromInfo.ReadToMap(chromFile)
	var lengthFromEnd int
	var found bool

	for i := range records {
		_, found = ref[records[i].Chrom]
		if found != true {
			log.Fatalf("Did not find '%s' in the chrom.sizes file", records[i].Chrom)
		}
		lengthFromEnd = ref[records[i].Chrom].Size - records[i].ChromEnd
		if lengthFromEnd < 0 {
			log.Fatalf("inputBed coordinates are outside chrom.sizes coordinate range, %s", records[i])
		}
		records[i].Score = numbers.Min(lengthFromEnd, records[i].ChromStart)
		if records[i].FieldsInitialized < 5 {
			records[i].FieldsInitialized = 5
		}
	}
	bed.Write(outFile, records)
}

func usage() {
	fmt.Print(
		"bedDistanceFromChrEnds - Returns a bed file with the Score field containing the minimum\n" +
			"distance from the end of the chromosome.\n" +
			"Usage:\n" +
			"bedDistanceFromEnds input.bed reference.chrom.sizes output.bed\n" +
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
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	chromFile := flag.Arg(1)
	outFile := flag.Arg(2)

	bedDistanceFromChrEnds(inFile, chromFile, outFile)
}
