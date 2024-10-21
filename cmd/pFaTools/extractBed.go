package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"log"
	"os"
)

// ExtractBedSettings defines the usage settings for the pFa extractBed subcommand.
type ExtractBedSettings struct {
	InFile     string
	Region     string
	OutFile    string
	TakeCoords bool
}

// extractBedUsage defines the usage statement for the pFaTools extractBed subcommand.
func extractBedUsage(extractBedFlags *flag.FlagSet) {
	fmt.Printf("pFaTools extractBed - a tool for extracting a subsequence from a pFa, as specified by a Bed file.\n" +
		"A note on input and output file types:\n" +
		"\tpFa (probabilistic FASTA) encodes position-wise probabilities in matrix entries. Columns should therefore sum to 1.\n" +
		"\tbed (browser extensible data) stores genomic regions as coordinates and associated annotations. User can specify multiple regions to be extracted." +
		"Usage:\n" +
		"pFaTools extractBed in.pFa region.bed outFile.pFa\n" +
		"options:\n")
	extractBedFlags.PrintDefaults()
}

// parseExtractBedArgs is the main function of the pFaTools extractBed subcommand. It parses options and runs the pFaExtractBed function.
func parseExtractBedArgs() {
	var expectedNumArgs int = 3
	var err error
	extractBedFlags := flag.NewFlagSet("extractBed", flag.ExitOnError)
	var takeCoords *bool = extractBedFlags.Bool("takeCoords", false, "Specify if bed coordinate fields should be used to differentiate regions in FASTA. Defaults to false.")
	err = extractBedFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	extractBedFlags.Usage = func() { extractBedUsage(extractBedFlags) }

	if len(extractBedFlags.Args()) != expectedNumArgs {
		extractBedFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(extractBedFlags.Args()))
	}

	inFile := extractBedFlags.Arg(0)
	region := extractBedFlags.Arg(1)
	outFile := extractBedFlags.Arg(2)

	s := ExtractBedSettings{
		InFile:     inFile,
		Region:     region,
		OutFile:    outFile,
		TakeCoords: *takeCoords,
	}

	pFaExtractBed(s)
}

// pFaExtractBed parses an input pFasta file and input Bed file, and extractBeds the file according to the Bed file.
func pFaExtractBed(s ExtractBedSettings) {
	records := pFasta.ExtractBed(pFasta.Read(s.InFile), bed.Read(s.Region), s.TakeCoords)
	pFasta.Write(s.OutFile, records)
}
