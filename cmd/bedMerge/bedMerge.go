// Command Group: "BED Tools"

// Combines overlapping bed entries, keeping max score. Output will be sorted by genome coordinate
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

func bedMerge(infile string, outfile string, mergeAdjacent int, lowMem bool, keepAllNames bool) {
	if lowMem {
		bedMergeLowMem(infile, outfile, mergeAdjacent)
	} else {
		bedMergeHighMem(infile, outfile, mergeAdjacent, keepAllNames)
	}
}

func bedMergeLowMem(infile string, outfile string, mergeAdjacent int) {
	var err error
	b := bed.GoReadToChan(infile)
	out := fileio.EasyCreate(outfile)
	var firstTime bool = true
	var currentMax bed.Bed
	var minDist int

	for i := range b {
		if firstTime {
			firstTime = false
			currentMax = i
		} else {
			minDist, err = bed.MinimumDistance(currentMax, i)
			if bed.Overlap(currentMax, i) || minDist <= mergeAdjacent && err == nil {
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

func bedMergeHighMem(infile string, outfile string, mergeAdjacent int, keepAllNames bool) {
	var records = bed.Read(infile)
	outList := bed.MergeHighMem(records, mergeAdjacent, keepAllNames)
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
	var pad *int = flag.Int("pad", -1, "Merge bed entries that are N bases away from each other. A pad value of 0 is equivalent to -mergeAdjacent")
	var lowMem *bool = flag.Bool("lowMem", false, "Use the low memory algorithm. Requires input file to be pre-sorted.")
	var keepAllNames *bool = flag.Bool("keepAllNames", false, "If set to true, merged beds will also have a merged name field in a comma separated list, cannot currently be combined with lowMem option")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	if *pad > -1 {
		*pad++
	}

	if *mergeAdjacent && *pad < 0 {
		*pad = 1
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)

	bedMerge(infile, outfile, *pad, *lowMem, *keepAllNames)
}
