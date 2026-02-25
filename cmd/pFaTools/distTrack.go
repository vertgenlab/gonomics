package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"math"
	"os"
	"strings"
)

// DistTrackSettings defines the usage settings for the pFa distTrack subcommand.
type DistTrackSettings struct {
	InFile1       string
	InFile2       string
	OutFile      string
	DefaultValue float64 // default value for wig.
	OutName  string
}

// distTrackUsage defines the usage statement for the pFaTools distTrack subcommand.
func distTrackUsage(distTrackFlags *flag.FlagSet) {
	fmt.Printf("pFaTools distTrack - a tool for plotting distance information between two sequences\n" +
		"uncertainty from an input pFa file.\n" +
		"Usage:\n" +
		"pFaTools distTrack A.pFa B.pFa out.wig" +
		"options:\n")
	distTrackFlags.PrintDefaults()
}

// parseDistTrackArgs is the main function of the pFaTools distTrack subcommand. It parses options and runs the pFaDistTrack function.
func parseDistTrackArgs() {
	var expectedNumArgs int = 3
	var err error
	distTrackFlags := flag.NewFlagSet("distTrack", flag.ExitOnError)
	var defaultValue *float64 = distTrackFlags.Float64("defaultValue", math.MaxFloat64, "Set the default value for the wig. Positions with this value will not be written to the file.")
	var outName *string = distTrackFlags.String("outName", "", "Specify name of the out file. If none is provided, output name will be the same as the first input name.")
	err = distTrackFlags.Parse(os.Args[3:])
	exception.PanicOnErr(err)
	distTrackFlags.Usage = func() { distTrackUsage(distTrackFlags) }
	if len(distTrackFlags.Args()) != expectedNumArgs {
		distTrackFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(distTrackFlags.Args()))
	}

	inFile1 := distTrackFlags.Arg(0)
	inFile2 := distTrackFlags.Arg(1)
	outFile := distTrackFlags.Arg(2)
	s := DistTrackSettings{
		InFile1:       inFile1,
		InFile2:       inFile2,
		OutFile:      outFile,
		DefaultValue: *defaultValue,
		OutName: *outName,
	}
	pFaDistTrack(s)
}

// pFaDistTrack reports a wig track two input singular pFastas (ie takes the first pfasta in the file), providing base-by-base
// information about the similarity of the pFastas. Assumes the pFastas are aligned.
func pFaDistTrack(s DistTrackSettings) {
	recordA := pFasta.Read(s.InFile1)
	var answer = make(map[string]wig.Wig, 1)
	if strings.HasSuffix(s.InFile2, ".fa") || strings.HasSuffix(s.InFile2, ".fasta") {
		recordB := fasta.Read(s.InFile2)
		answer[s.OutName] = pFasta.DistTrackFasta(recordA[0], recordB[0], s.OutName, s.DefaultValue)
	} else {
		recordB := pFasta.Read(s.InFile2)
		answer[s.OutName] = pFasta.DistTrack(recordA[0], recordB[0], s.OutName, s.DefaultValue)
	}

	
	
	wig.Write(s.OutFile, answer)
}
