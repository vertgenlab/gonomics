package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"log"
)

func bedpeOverlap(bedpeSelectFile string, bedpeInFile string, contactOutFile string) {
	var inIntervals = make([]interval.Interval, 0)
	var leftOverlaps []interval.Interval
	var rightOverlaps []interval.Interval
	var rightHalf bedpe.BedPeHalf
	var leftHalf bedpe.BedPeHalf
	var err error

	out := fileio.EasyCreate(contactOutFile)

	contactRecords := bedpe.Read(bedpeSelectFile)

	inBedPe := bedpe.Read(bedpeInFile)

	for _, i := range inBedPe {
		left, right := bedpe.SplitBedPe(i)
		inIntervals = append(inIntervals, left, right)
	}
	inTree := interval.BuildTree(inIntervals)

	for _, i := range contactRecords {
		leftOverlaps = interval.Query(inTree, i.A, "any")
		rightOverlaps = interval.Query(inTree, i.B, "any")
		found := false
		for j := range leftOverlaps {
			for k := range rightOverlaps {
				leftHalf = leftOverlaps[j].(bedpe.BedPeHalf)
				rightHalf = rightOverlaps[k].(bedpe.BedPeHalf)
				if leftHalf.Home == rightHalf.Home {
					found = true
				}
			}
		}
		if found {
			bedpe.WriteToFileHandle(out, *leftHalf.Home)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)

}

func usage() {
	fmt.Print("bedpeOverlap - Returns all bedpe entries that fully overlap another bedpe entry. " +
		"bedpe entries from in the inFile will be selected based on overlaps with the selectFile.\n" +
		"Usage:\n" +
		"	bedpeOverlap [options] selectFile.bedpe inFile.bedPe out.bedpe\n\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
	}

	bedpeSelectFile := flag.Arg(0)
	bedpeInFile := flag.Arg(1)
	contactOutFile := flag.Arg(2)

	bedpeOverlap(bedpeSelectFile, bedpeInFile, contactOutFile)
}
