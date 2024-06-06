package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"math"
	"os"
)

// EntropyTrackSettings defines the usage settings for the pFa entropyTrack subcommand.
type EntropyTrackSettings struct {
	InFile       string
	OutFile      string
	DefaultValue float64 // default value for wig.
}

// entropyTrackUsage defines the usage statement for the pFaTools entropyTrack subcommand.
func entropyTrackUsage(entropyTrackFlags *flag.FlagSet) {
	fmt.Printf("pFaTools entropyTrack - a tool for plotting positional information on sequence\n" +
		"uncertainty from an input pFa file.\n" +
		"Usage:\n" +
		"pFaTools entropyTrack in.pFa out.wig" +
		"options:\n")
	entropyTrackFlags.PrintDefaults()
}

// parseEntropyTrackArgs is the main function of the pFaTools entropyTrack subcommand. It parses options and runs the pFaEntropyTrack function.
func parseEntropyTrackArgs() {
	var expectedNumArgs int = 2
	var err error
	entropyTrackFlags := flag.NewFlagSet("entropyTrack", flag.ExitOnError)
	var defaultValue *float64 = entropyTrackFlags.Float64("defaultValue", math.MaxFloat64, "Set the default value for the wig. Positions with this value will not be written to the file.")
	err = entropyTrackFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	entropyTrackFlags.Usage = func() { entropyTrackUsage(entropyTrackFlags) }
	if len(entropyTrackFlags.Args()) != expectedNumArgs {
		entropyTrackFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(entropyTrackFlags.Args()))
	}

	inFile := entropyTrackFlags.Arg(0)
	outFile := entropyTrackFlags.Arg(1)
	s := EntropyTrackSettings{
		InFile:       inFile,
		OutFile:      outFile,
		DefaultValue: *defaultValue,
	}
	pFaEntropyTrack(s)
}

// pFaEntropyTrack reports a wig track from an input pFasta file, providing base-by-base
// information about the shannon entropy of the pBase, reflecting the sequence uncertainty.
func pFaEntropyTrack(s EntropyTrackSettings) {
	var currPos int
	var currentWig wig.Wig
	records := pFasta.Read(s.InFile)
	var answer = make(map[string]wig.Wig, len(records))

	for _, v := range records {
		currentWig = wig.Wig{StepType: "fixedStep", Chrom: v.Name, Start: 1, Step: 1, DefaultValue: s.DefaultValue}
		currentWig.Values = make([]float64, len(v.Seq))
		for currPos = 0; currPos < len(v.Seq); currPos++ {
			currentWig.Values[currPos] = pDna.Entropy(v.Seq[currPos])
		}
		answer[v.Name] = currentWig
	}
	wig.Write(s.OutFile, answer)
}
