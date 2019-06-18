package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"fmt"
	"log"
	"github.com/vertgenlab/gonomics/common"
	"flag"
)

func bedMerge(infile string, outfile string) {
	var records []*bed.Bed = bed.Read(infile)
	var outlist []*bed.Bed 
	var currentMax *bed.Bed = records[0]

	//assumes sorted bed


	for i := 0; i < len(records); i++ {
		if overlap(currentMax, records[i]) {
			if records[i].Score > currentMax.Score {
				currentMax.Score = records[i].Score
			}
			currentMax.ChromEnd = common.MaxInt64(records[i].ChromEnd, currentMax.ChromEnd)
		} else {
			outlist = append(outlist, currentMax)
			currentMax = records[i]
		}
	}
	outlist = append(outlist, currentMax)
	bed.Write(outfile, outlist, 5)
}

func overlap(bed1 *bed.Bed, bed2 *bed.Bed) bool {
	if bed1.Chrom != bed2.Chrom {
		return false
	} else if bed1.ChromEnd < bed2.ChromStart || bed2.ChromEnd < bed1.ChromStart {
		return false
	}
	return true
}

func usage() {
	fmt.Print(
                "bedMerge - Combines overlapping bed entries, keeping max score.\n" +
                        "Usage:\n" +
                        "bedMerge input.bed output.bed\n" +
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

	infile := flag.Arg(0)
	outfile := flag.Arg(1)

	bedMerge(infile, outfile)
}