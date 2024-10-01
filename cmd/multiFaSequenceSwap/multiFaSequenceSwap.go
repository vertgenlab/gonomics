package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
)

// multiFaSubsequenceSwap swaps a series of specified regions from a bed file between two sequences in a MultiFa object
// The foregroundName is the name of the sequence containing the regions in the bed elements;  backgroundName is the name of the sequence including the regions outside the bed elements
// input of this function includes the multiFa file name, the index of the background and foreground sequences from the multiFa, and a bed object with the relevant
// start and stop indices of the swap regions
// it will return a file containing a Fasta slice (which includes several sequences), with one more sequence that the sequences in the original
// multiFa; this extra sequence is the background sequence with relevant foreground spots swapped as dictated by the bed file.
func multiFaSubsequenceSwap(inFile string, swapRegionsFile string, backgroundName, foregroundName string, chromName string, outFile string) {
	var currRefPos, currAlnPos, lastRefPos, lastAlnPos int //these variables will be used to swap regions later in the code
	// Load the original sequences from the multiFa file.
	records := fasta.Read(inFile)
	swapRegions := bed.Read(swapRegionsFile)
	bed.SortByCoord(swapRegions)
	swapRegions = bed.MergeBeds(swapRegions) // avoids the case of overlapping foreground regions
	background := getFaIndex(records, backgroundName)
	foreground := getFaIndex(records, foregroundName)
	// Ensure the specified background and foreground sequence indices are within the valid range.
	if background < 0 || background >= len(records) || foreground < 0 || foreground >= len(records) {
		log.Fatalf("Error: Invalid sequence index.\n")

	}
	answerSeq := fasta.Copy(records[background]) //copy of the background sequence is made
	for _, currRegion := range swapRegions {     //loop over every single swap region in the bed
		// Check if swap region is on the provided chromsome, otherwise skip iteration
		if currRegion.Chrom != chromName {
			continue
		}

		// Validate the swap region.
		if currRegion.ChromStart < 0 || currRegion.ChromStart >= currRegion.ChromEnd {
			log.Fatalf("Error: Invalid swap region. \n")
		}

		// we've merged and sorted beds by position, so we should always start below the chromStart in ref positions.
		if currRefPos > currRegion.ChromStart {
			log.Fatalf("Something went wrong. Debug required.\n")
		}

		for currRefPos < currRegion.ChromEnd {
			if currRefPos >= currRegion.ChromStart {
				answerSeq.Seq[currAlnPos] = records[foreground].Seq[currAlnPos]
			}
			currAlnPos++
			currRefPos = fasta.AlnPosToRefPosCounter(records[0], currAlnPos, lastRefPos, lastAlnPos)
			lastRefPos, lastAlnPos = currRefPos, currAlnPos
		}
	}
	var finalSeqs []fasta.Fasta
	answerSeq.Name = fmt.Sprintf("%v.swapped", backgroundName)
	finalSeqs = append(records, answerSeq)
	fasta.Write(outFile, finalSeqs)

}
func getFaIndex(records []fasta.Fasta, name string) int {
	for currIndex := 0; currIndex < len(records); currIndex++ {
		if records[currIndex].Name == name {
			return currIndex
		}
	}
	log.Fatalf("Error: The requested sequence name: %s does not exist in the provided input. \n", name)
	return -1

}

func usage() {
	fmt.Print(
		"multiFaSequenceSwap - Swaps regions of a sequence with regions from a different corresponding sequence. \n" +
			"Input bed regions should reflect reference coordinates for the multiFa.\n" +
			"foregroundName is the name of the sequence that includes the regions in the bed elements.\n" +
			"backgroundName is the name of the sequence that includes the regions outside the bed elements.\n" +
			"Input of this function includes the multiFa file name,\n" +
			"the index of the background and foreground sequences from the multiFa,\n" +
			"a bed object with the relevant start and stop indices of the swap regions,\n" +
			"the name of the outfile multiFa, and the name of the resulting sequence.\n" +
			"It will return a file containing one more multiFa than is in the inFile,\n" +
			"named backgroundName.swapped.\n" +
			"Usage:  \n" +
			"multiFaSequenceSwap inFile.fa bedFile.bed background foreground chromName outFile.fa \n" +
			"options: \n",
	)
	flag.PrintDefaults()

}

func main() {

	var expectedNumArgs int = 6
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	bedFile := flag.Arg(1)
	backgroundName := flag.Arg(2)
	foregroundName := flag.Arg(3)
	chromName := flag.Arg(4)
	outFile := flag.Arg(5)

	multiFaSubsequenceSwap(inFile, bedFile, backgroundName, foregroundName, chromName, outFile)

}
