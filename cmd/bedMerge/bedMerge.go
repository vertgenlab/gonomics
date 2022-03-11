// Command Group: "BED Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

func bedMerge(infile string, outfile string, mergeAdjacent bool, lowMem bool) {
	if lowMem {
		bedMergeLowMem(infile, outfile, mergeAdjacent)
	} else {
		bedMergeHighMem(infile, outfile, mergeAdjacent)
	}
}

func bedMergeLowMem(infile string, outfile string, mergeAdjacent bool) {
	var err error
	b := bed.GoReadToChan(infile)
	out := fileio.EasyCreate(outfile)
	var firstTime bool = true
	var currentMax bed.Bed

	for i := range b {
		if firstTime {
			firstTime = false
			currentMax = i
		} else {
			if bed.Overlap(currentMax, i) || mergeAdjacent && bed.Adjacent(currentMax, i) {
				if i.Score > currentMax.Score {
					currentMax.Score = i.Score
				}
				currentMax.ChromEnd = numbers.Max(i.ChromEnd, currentMax.ChromEnd)
			} else {
				bed.WriteBed(out, currentMax)
				currentMax = i
			}
		}
	}
	bed.WriteBed(out, currentMax)

	err = out.Close()
	exception.PanicOnErr(err)
}

func bedMergeHighMem(infile string, outfile string, mergeAdjacent bool) {
	var records []bed.Bed = bed.Read(infile)
	var outList []bed.Bed
	var currentMax bed.Bed = records[0]

	bed.SortByCoord(records)

	for i := 1; i < len(records); i++ {
		if bed.Overlap(currentMax, records[i]) || mergeAdjacent && bed.Adjacent(currentMax, records[i]) {
			if records[i].Score > currentMax.Score {
				currentMax.Score = records[i].Score
			}
			currentMax.ChromEnd = numbers.Max(records[i].ChromEnd, currentMax.ChromEnd)
		} else {
			outList = append(outList, currentMax)
			currentMax = records[i]
		}
	}
	outList = append(outList, currentMax)
	bed.Write(outfile, outList)
}

func usage() {
	fmt.Print(
		"bedMerge - Combines overlapping bed entries, keeping max score. Output will be sorted by genome coordinate.\n" +
			"Usage:\n" +
			"bedMerge input.bed output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var mergeAdjacent *bool = flag.Bool("mergeAdjacent", false, "Merge non-overlapping entries with direct adjacency.")
	var lowMem *bool = flag.Bool("lowMem", false, "Use the low memory algorithm. Requires input file to be pre-sorted.")

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

	bedMerge(infile, outfile, *mergeAdjacent, *lowMem)
}
