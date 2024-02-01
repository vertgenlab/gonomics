package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

// multiFaSubsequenceSwap swaps a series of specified regions from a bed file between two sequences in a MultiFa object
// The foregroundName is the name of the sequence containing the regions in the bed elements;  backgroundName is the name of the sequence including the regions outside the bed elements
// input of this function includes the multiFa file name, the index of the background and foreground sequences from the multiFa, and a bed object with the relevant
// start and stop indices of the swap regions
// it will return a file containing a Fasta slice (which includes several sequences), with one more sequence that the sequences in the original
// multiFa; this extra sequence is the background sequence with relevant foreground spots swapped as dictated by the bed file.
func multiFaSubsequenceSwap(inFile string, swapRegionsFile string, backgroundName, foregroundName string, outFile string, outSeqName string) {
	var currStart, currEnd, currPos int //these variables will be used to swap regions later in the code
	// Load the original sequences from the multiFa file.
	originalSeqs := fasta.Read(inFile)
	swapRegions := bed.Read(swapRegionsFile)
	background := getFaIndex(originalSeqs, backgroundName)
	foreground := getFaIndex(originalSeqs, foregroundName)
	// Ensure the specified background and foreground sequence indices are within the valid range.
	if background < 0 || background >= len(originalSeqs) || foreground < 0 || foreground >= len(originalSeqs) {
		log.Fatalf("Error: Invalid sequence index.\n")

	}
	answerSeq := fasta.Copy(originalSeqs[background]) //copy of the background sequence is made
	for _, swapRegion := range swapRegions {          //loop over every single swap region in the bed
		currStart = swapRegion.ChromStart
		currEnd = swapRegion.ChromEnd

		// Validate the swap region.
		if currStart < 0 || currEnd > len(answerSeq.Seq) || currStart >= currEnd {
			log.Fatalf("Error: Invalid swap region. \n")

		}

		// Swap the corresponding portion of answerSeq with the foreground sequence.
		for currPos = currStart; currPos < currEnd; currPos++ {
			answerSeq.Seq[currPos] = originalSeqs[foreground].Seq[currPos]
		}
	}
	var finalSeqs []fasta.Fasta
	answerSeq.Name = outSeqName
	finalSeqs = append(originalSeqs, answerSeq)

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
			"foregroundName is the name of the sequence that includes the regions in the bed elements.\n" +
			"backgroundName is the name of the sequence that includes the regions outside the bed elements.\n" +
			"Input of this function includes the multiFa file name,\n" +
			"the index of the background and foreground sequences from the multiFa,\n" +
			"a bed object with the relevant start and stop indices of the swap regions,\n" +
			"the name of the outfile multiFa, and the name of the resulting sequence.\n" +
			"It will return a file containing one more multiFa than is in the inFile.\n" +
			"Usage:  \n" +
			"multiFaSequenceSwap inFile.fa bedFile.bed background foreground outFile.fa outSeqName \n" +
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
	outFile := flag.Arg(4)
	outSeqName := flag.Arg(5)

	multiFaSubsequenceSwap(inFile, bedFile, backgroundName, foregroundName, outFile, outSeqName)

}
