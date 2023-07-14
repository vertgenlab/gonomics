// Command Group: "BED Tools"

// Returns a bed file with the Score field containing the minimum distance from the end of the chromosome.
package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/numbers"
)

func bedDistanceFromChrEnds(inFile string, chromFile string, outFile string) {
	records := bed.Read(inFile)           //read the bed infile and store it in records
	ref := chromInfo.ReadToMap(chromFile) //access the information in chrom.sizes file for the genome the bed originates from
	var lengthFromEnd int
	var found bool

	for i := range records { //for each bed
		_, found = ref[records[i].Chrom] // check if the chrom in the bed exists in the chrom.sizes file
		if !found {                      //if the chrom doesn't exist, throw an error
			log.Fatalf("Did not find '%s' in the chrom.sizes file", records[i].Chrom)
		}
		lengthFromEnd = ref[records[i].Chrom].Size - records[i].ChromEnd //calculates the distance from the right termination of the chrom to the end of this bed record
		if lengthFromEnd < 0 {                                           //throw an error if the length is less than 0
			log.Fatalf("inputBed coordinates are outside chrom.sizes coordinate range, %s", records[i])
		}
		records[i].Score = numbers.Min(lengthFromEnd, records[i].ChromStart) //calculate if the beginning of the bed is closer to the start than the end of the bed is to the end of the chrom and store that distance in the score field of the bed that will be written out
		if records[i].FieldsInitialized < 5 {                                // output bed should have 5 fields (more information about the bed fields in the codebase can be found in gonomics/bed/bed.go
			records[i].FieldsInitialized = 5
		}
	}
	bed.Write(outFile, records) //write the beds with their Score field containing their distance to the end of the chromosome out to the specified file
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
