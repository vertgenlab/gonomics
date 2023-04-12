package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/interval/lift"
	"log"
)

func bedpeOverlap(selectFile string, bedpeInFile string, contactOutFile string, bedSelect bool, overlapThreshold float64, overlapBoth bool) {
	//bedpeOverlap will work with either a bedpe select file or a bed select file. First we determine which program to run.
	if bedSelect && overlapBoth {
		SelectIsBedBoth(selectFile, bedpeInFile, overlapThreshold, contactOutFile)
	} else if bedSelect {
		SelectIsBed(selectFile, bedpeInFile, overlapThreshold, contactOutFile)
	} else {
		SelectIsBedPe(selectFile, bedpeInFile, contactOutFile)
	}
}

// function to calculate what percent of a bedpeHalf overlaps the selectBed
func overlapPercent(possOverlaps interval.Interval, halfBedPe bed.Bed) float64 {
	var overlapEntryStart, overlapEntryEnd, halfBedPeStart, halfBedPeEnd, overlapSize int
	var answer float64

	overlapEntryStart = possOverlaps.GetChromStart()
	overlapEntryEnd = possOverlaps.GetChromEnd()
	halfBedPeStart = halfBedPe.ChromStart
	halfBedPeEnd = halfBedPe.ChromEnd

	overlapSize = lift.MatchOverlapLen(overlapEntryStart, overlapEntryEnd, halfBedPeStart, halfBedPeEnd)

	answer = float64(overlapSize) / float64(halfBedPeEnd-halfBedPeStart)
	return answer
}

///func overlapPercentBedPeSelect (leftSideA interval.Interval, leftSideB interval.Interval, rightSideA interval.Interval, rightSideB interval.Interval) float64 {
//	var rightSideStartA int
//	var rightSideEndA int
//	var rightSideStartB int
//	var rightSideEndB int
//	var leftSideStartA int
//	var leftSideEndA int
//	var leftSideStartB int
//	var leftSideEndB int
//	var overlapSizeRight int
//	var overlapSizeLeft int
//	var overlapPercentLeft float64
//	var overlapPercentRight float64
//
//	rightSideStartA = rightSideA.GetChromStart()
//	rightSideEndA = rightSideA.GetChromEnd()
//	rightSideStartB = rightSideB.GetChromStart()
//	rightSideEndB = rightSideB.GetChromEnd()
//
//	leftSideStartA = leftSideA.GetChromStart()
//	leftSideEndA = leftSideA.GetChromEnd()
//	leftSideStartB = leftSideB.GetChromStart()
//	leftSideEndB = leftSideB.GetChromEnd()
//
//	overlapSizeRight = lift.MatchOverlapLen(rightSideStartA, rightSideEndA, rightSideStartB, rightSideEndB)
//
//	overlapSizeLeft = lift.MatchOverlapLen(leftSideStartA, leftSideEndA, leftSideStartB, leftSideEndB)
//
//need to resolve this by integrating which is the select bedpe vs the in bedpe. in bedpe should be the denominator?
//	overlapPercentRight = float64(overlapSizeRight) / float64()
//	overlapPercentLeft = float64(overlapSizeRight) / float64()
//
//	return math.Min(overlapPercentRight, overlapPercentLeft)
//}
///

// This function is for the case where the select file is a bed.
// input bedpe entries are retained if either end overlaps one of the bedSelectFile entries.
func SelectIsBed(bedSelectFile string, bedpeInFile string, overlapThreshold float64, contactOutFile string) {
	var selectIntervals = make([]interval.Interval, 0)
	var currOverlaps []interval.Interval
	var err error
	var found bool

	selectRecords := bed.Read(bedSelectFile)
	inBedPe := bedpe.Read(bedpeInFile)
	out := fileio.EasyCreate(contactOutFile)

	for _, i := range selectRecords {
		selectIntervals = append(selectIntervals, i)
	}
	selectTree := interval.BuildTree(selectIntervals)

	for _, i := range inBedPe {
		currOverlaps = interval.Query(selectTree, i.A, "any")
		if len(currOverlaps) > 0 { //if A, the left side of the input bedpe, overlaps any of the select beds, write the bedpe to output.
			if overlapThreshold == 0 {
				bedpe.WriteToFileHandle(out, i)
			} else {
				found = false
				for _, j := range currOverlaps {
					if !found && overlapPercent(j, i.A) >= overlapThreshold {
						found = true
						bedpe.WriteToFileHandle(out, i)
					}
				}
			}
		} else { //otherwise check the right side (B)
			currOverlaps = interval.Query(selectTree, i.B, "any")
			if len(currOverlaps) > 0 {
				if overlapThreshold == 0 {
					bedpe.WriteToFileHandle(out, i)
				} else {
					found = false
					for _, j := range currOverlaps {
						if !found && overlapPercent(j, i.B) >= overlapThreshold {
							bedpe.WriteToFileHandle(out, i)
						}
					}
				}
			}
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func SelectIsBedBoth(bedSelectFile string, bedpeInFile string, overlapThreshold float64, contactOutFile string) {
	var selectIntervals = make([]interval.Interval, 0)
	var aOverlaps []interval.Interval
	var bOverlaps []interval.Interval
	var err error
	var found bool

	selectRecords := bed.Read(bedSelectFile)
	inBedPe := bedpe.Read(bedpeInFile)
	out := fileio.EasyCreate(contactOutFile)

	for _, i := range selectRecords {
		selectIntervals = append(selectIntervals, i)
	}
	selectTree := interval.BuildTree(selectIntervals)

	var j, k interval.Interval
	var i bedpe.BedPe

	for _, i = range inBedPe {
		aOverlaps = interval.Query(selectTree, i.A, "any")
		if len(aOverlaps) > 0 { //if A, the left side of the input bedpe, overlaps any of the select beds, write the bedpe to output.
			if overlapThreshold == 0 {
				bOverlaps = interval.Query(selectTree, i.B, "any")
				if len(bOverlaps) > 0 {
					bedpe.WriteToFileHandle(out, i)
				}
			} else {
				found = false
				for _, j = range aOverlaps {
					if !found && overlapPercent(j, i.A) >= overlapThreshold {
						bOverlaps = interval.Query(selectTree, i.B, "any")
						for _, k = range bOverlaps {
							if !found && overlapPercent(k, i.B) >= overlapThreshold {
								found = true
								bedpe.WriteToFileHandle(out, i)
							}
						}
					}
				}
			}
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// This is the case where the select file is a bedpe. Input bedpe entries will be retained
// in the output if both ends of the bedpe overlap both ends of a bedpe entry in the select file.
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
		"When the select file is a bed, as specified in the option 'bedSelect',\n" +
		"entries are retained if at least one end of an input bedpe overlaps a bed entry\n" +
		"in the select file.\n" +
		"overlapThreshold will only return entries that have an overlap that covers greater than x% of the bedpe overlap entry\n" +
		"by default every overlap regardless of overlap percentage will be reported\n" +
		"overlapThreshold is only compatible with -bedSelect\n" +
		"Usage:\n" +
		"bedpeOverlap [options] selectFile inputFile.bedpe out.bedpe\n\n")
	flag.PrintDefaults()
}

func main() {
	var bedSelect *bool = flag.Bool("bedSelect", false, "Set select file to be a BED file instead of a bedpe.")
	var overlapThreshold *float64 = flag.Float64("overlapThreshold", 0, "threshold that the percent overlap of a bedpe half to the select bed must satisfy. Must be a value between 0 and 1.")
	var overlapBoth *bool = flag.Bool("overlapBoth", false, "restricts outputs of -bedSelect to bedpe entries where both ends overlap a selectBed entry")

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

	if *overlapThreshold < 0 || *overlapThreshold > float64(1) {
		log.Fatalf("Error: overlap threshold must be between 0 and 1")
	}

	//should no longer be needed when overlapThresholdBedPe is implemented
	if *overlapThreshold != 0 && *bedSelect == false {
		log.Fatalf("Error: overlapThreshold must be used with bedSelect")
	}

	if *overlapBoth == true && *bedSelect == false {
		log.Fatalf("Error: overlapBoth must be used with bedSelect")
	}

	bedpeOverlap(SelectFile, bedpeInFile, contactOutFile, *bedSelect, *overlapThreshold, *overlapBoth)
}
