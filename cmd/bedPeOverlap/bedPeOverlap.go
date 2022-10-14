package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"log"
)

func bedpeOverlap(selectFile string, bedpeInFile string, contactOutFile string, bedSelect bool) {
	if bedSelect {
		SelectIsBed(selectFile, bedpeInFile, contactOutFile)
	} else {
		SelectIsBedPe(selectFile, bedpeInFile, contactOutFile)
	}
}

func SelectIsBed(bedSelectFile string, bedpeInFile string, contactOutFile string) {
	var selectIntervals = make([]interval.Interval, 0)
	var currOverlaps []interval.Interval
	var err error

	selectRecords := bed.Read(bedSelectFile)
	inBedPe := bedpe.Read(bedpeInFile)
	out := fileio.EasyCreate(contactOutFile)

	for _, i := range selectRecords {
		selectIntervals = append(selectIntervals, i)
	}
	selectTree := interval.BuildTree(selectIntervals)

	for _, i := range inBedPe {
		currOverlaps = interval.Query(selectTree, i.A, "any")
		if len(currOverlaps) > 0 {//if A, the left side of the input bedPe, overlaps any of the select beds, write the bedPe to output.
			bedpe.WriteToFileHandle(out, i)
		} else { //otherwise check the right side (B)
			currOverlaps = interval.Query(selectTree, i.B, "any")
			if len(currOverlaps) > 0 {
				bedpe.WriteToFileHandle(out, i)
			}
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func SelectIsBedPe(bedpeSelectFile string, bedpeInFile string, contactOutFile string) {
	var inIntervals = make([]interval.Interval, 0)
	var leftOverlaps, rightOverlaps []interval.Interval
	var rightHalf, leftHalf bedpe.BedPeHalf
	var err error

	out := fileio.EasyCreate(contactOutFile)

	contactRecords := bedpe.Read(bedpeSelectFile)

	inBedPe := bedpe.Read(bedpeInFile)

	var left, right bedpe.BedPeHalf
	for _, i := range inBedPe {
		left, right = bedpe.SplitBedPe(i)
		inIntervals = append(inIntervals, left, right)
	}
	inTree := interval.BuildTree(inIntervals)

	var found bool
	var j, k int
	for _, i := range contactRecords {
		found = false
		leftOverlaps = interval.Query(inTree, i.A, "any")
		rightOverlaps = interval.Query(inTree, i.B, "any")
		for j = range leftOverlaps {
			for k = range rightOverlaps {
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
	fmt.Print("bedpeOverlap - Filters bedpe entries based on overlaps from the select file.\n" +
		"Default behavior expects a bedpe select file and returns entries where both ends of a bedpe entry from the input file" +
		"overlap both ends of a bedpe entry from the select file.\n" +
		"Usage:\n" +
		"bedpeOverlap [options] selectFile inputFile.bedPe out.bedpe\n\n")
	flag.PrintDefaults()
}

func main() {
	var bedSelect *bool = flag.Bool("bedInput", false, "Set select file to be a BED file instead of a bedpe.")

	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
	}

	SelectFile := flag.Arg(0)
	bedpeInFile := flag.Arg(1)
	contactOutFile := flag.Arg(2)

	bedpeOverlap(SelectFile, bedpeInFile, contactOutFile, *bedSelect)
}
