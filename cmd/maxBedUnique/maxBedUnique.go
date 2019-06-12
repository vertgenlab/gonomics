package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"fmt"
	"log"
	"flag"
)

func maxBedUnique (infile string, outfile string) {
	var records []*bed.Bed = bed.Read(infile)
	var outlist []*bed.Bed
	var currentMax *bed.Bed = records[0]
	var goingUp bool = true
	var currentMin *bed.Bed = records[0]

	for i := 0; i < len(records); i++ {
		if goingUp {
			if bed.Overlap(currentMax, records[i]) {
				if (currentMax.Score < records[i].Score) {
                                currentMax = records[i]
				}
			} else {
				outlist = append(outlist, currentMax)
				currentMax = records[i]
				goingUp = false
				currentMin = currentMax
			}
		} else if (records[i].Score > currentMin.Score) {
			goingUp = true
			currentMax = records[i]
		}
	}

	outlist = append(outlist, currentMax)
	bed.Write(outfile, outlist, 5) //third input species field number
}

func usage() {
        fmt.Print(
                "maxBedUnique\n" +
                        "Usage:\n" +
                        "maxBedUnique input.bed output.bed\n" +
                        "options:\n")
        flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

        if len(flag.Args()) != expectedNumArgs {
                flag.Usage()
                log.Fatalf("Error: expecting %d arguments, but got %d\n",
                        expectedNumArgs, len(flag.Args()))
        }

        inFile := flag.Arg(0)
        outFile := flag.Arg(1)

	maxBedUnique(inFile, outFile)
}
