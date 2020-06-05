package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/sort"
	"log"
)

func usage() {
	fmt.Print(
		"mergesort - Executes an external merge sort of the input file based on desired sort criteria. \n" +
			"\t The input file should have a proper file extension depending on the input file type." +
			"\n\t Compatible input types currently include AXT, BED, VCF, SAM, and GIRAF\n" +
			"Usage:\n" +
			" mergesort [options] input.filetype outputFile\n")
	flag.PrintDefaults()
}

func mergeSort(inFile string, outFile string, tmpFilePrefix string, numLinesPerChunk int, sortCriteria string) {
	sort.ExternalMergeSort(inFile, numLinesPerChunk, tmpFilePrefix, outFile, sortCriteria)
}

func main() {
	expectedNumArgs := 2
	var tmpFilePrefix *string = flag.String("prefix", "tmp", "Prefix to use when writing temporary files (e.g. tmp_0).``")
	var numLinesPerChunk *int = flag.Int("tmpsize", 1000000, "The number of records to read into memory before writing to a tmp file.``")
	var byGenomicCoordinates *bool = flag.Bool("genomicCoordinates", true, "Sort Criteria: Chromosome -> StartPos -> EndPos.``")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("ERROR: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	// Default sort order is by genomic coordinates
	var sortCriteria string = "byGenomicCoordinates"

	var numSortCriteria int = 0
	// if newSortCriteria {
	// 		numSortCriteria++
	// }                         // ... for all sort criteria

	if numSortCriteria > 1 {
		log.Fatalln("Must use exactly 1 sort criteria.")
	}

	switch {
	// case newSortCriteria:
	//		sortCriteria = "newSortCriteria"
	case *byGenomicCoordinates: // Always last case
		sortCriteria = "byGenomicCoordinates"
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	mergeSort(inFile, outFile, *tmpFilePrefix, *numLinesPerChunk, sortCriteria)
}
