package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
)

// SubsequenceSwap swaps a series of specified regions from a bed file between two sequences in a MultiFa object
// The sequence called "background" is taken, and parts of "foreground" are swapped with the corresponding parts in background
// input of this function includes the multiFa file name, the index of the background and foreground sequences from the multiFa, and a bed object with the relevant
// start and stop indices of the swap regions
// it will return a file containing a Fasta slice (which includes several sequences), with one more sequence that the sequences in the original
// multiFa; this extra sequence is the background sequence with relevant foreground spots swapped as dictated by the bed file.
func SubsequenceSwap(inFile string, swapRegionsFile string, background, foreground int, outFile string, outSeqName string) {
	var currStart, currEnd, currPos int //these variables will be used to swap regions later in the code

	// Load the original sequences from the multiFa file.
	originalSeqs := fasta.Read(inFile) //Converting a multiFa to a Fasta file
	lines := fileio.Read(swapRegionsFile)
	fmt.Println(lines)
	swapRegions := bed.Read(swapRegionsFile)

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

	// Create and return the new Fasta with one more sequence than the input.
}
func usage() {
	fmt.Print(
		"multiFaSequenceSwap - Swaps regions of a sequence with regions from a different corresponding sequence. \n" +
			"Usage: \n" +
			"options: \n",
	)
	flag.PrintDefaults()

}

func main() {

	var expectedNumArgs int = 5
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d/n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	bedFile := flag.Arg(1)
	background := parse.StringToInt(flag.Arg(2))
	foreground := parse.StringToInt(flag.Arg(3))
	outFile := flag.Arg(4)
	outSeqName := flag.Arg(5)

	SubsequenceSwap(inFile, bedFile, background, foreground, outFile, outSeqName)

}
