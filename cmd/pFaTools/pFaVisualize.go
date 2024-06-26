package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/browser"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"os"
	"strings"
)

// VisualizeSettings defines the usage settings for the pFa Visualize subcommand.
type VisualizeSettings struct {
	InFile           string
	OutDir           string
	Start            int
	End              int
	SigFigs          int
	DecimalPlaces    int
	LineLength       int
	Chrom            string
	StartOfAlignment bool
	EndOfAlignment   bool
	TsvOut           bool
}

// VisualizeUsage defines the usage statement for the pFaTools Visualize subcommand.
func VisualizeUsage(VisualizeFlags *flag.FlagSet) {
	fmt.Printf("pFaTools Visualize - Provides human-readable sequence from a given pFa.\n" +
		"Start and end arguments are 0-indexed. The visualiser is inclusive of the start \n" +
		"index, exclusive of the end index.\n" +
		"Keyword 'START' for the start argument makes a visualization from the start of the pfasta.\n" +
		"Keyword 'END' for the end argument makes a visualization until the end of the pfasta.\n" +
		"Usage:\n" +
		"PFaTools visualize in.pfa start end outDir.txt\n" +
		"options:\n")
	VisualizeFlags.PrintDefaults()
}

// parseVisualizeArgs is the main function of the pFaTools Visualize subcommand. It parses options and runs the pFaVisualize function.
func parseVisualizeArgs() {
	var expectedNumArgs int = 4
	var err error
	VisualizeFlags := flag.NewFlagSet("Visualize", flag.ExitOnError)
	var sigFigs *int = VisualizeFlags.Int("sigFigs", 0, "Specify the number of significant figures to round to. Leave as 0 if decimal rounding desired.")
	var decimalPlaces *int = VisualizeFlags.Int("decimal", 5, "Specify the number of decimal places to round to. Should be between 1 and 15. Note that decimal places more to the right will accrue larger floating point error.")
	var lineLength *int = VisualizeFlags.Int("lineLength", 50, "Sets length of each alignment line.")
	var chrom *string = VisualizeFlags.String("chrom", "", "Specify the name of the sequence to display. Can be empty if only one sequence in input pfasta.")
	var tsvOut *bool = VisualizeFlags.Bool("tsvOut", false, "Specify if user wants R-ready output (true). Defaults to false.")
	var startOfAlignment bool = false
	var endOfAlignment bool = false
	var start int
	var end int

	err = VisualizeFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	VisualizeFlags.Usage = func() { VisualizeUsage(VisualizeFlags) }

	if len(VisualizeFlags.Args()) != expectedNumArgs {
		VisualizeFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(VisualizeFlags.Args()))
	}

	inFile := VisualizeFlags.Arg(0)
	if strings.ToLower(VisualizeFlags.Arg(1)) == "start" {
		startOfAlignment = true
		start = 0
	} else {
		start = parse.StringToInt(VisualizeFlags.Arg(1))
	}
	if strings.ToLower(VisualizeFlags.Arg(2)) == "end" {
		endOfAlignment = true
		end = -1
	} else {
		end = parse.StringToInt(VisualizeFlags.Arg(2))
	}
	outDir := VisualizeFlags.Arg(3)

	s := VisualizeSettings{
		InFile:           inFile,
		OutDir:           outDir,
		Start:            start,
		End:              end,
		SigFigs:          *sigFigs,
		DecimalPlaces:    *decimalPlaces,
		LineLength:       *lineLength,
		Chrom:            *chrom,
		StartOfAlignment: startOfAlignment,
		EndOfAlignment:   endOfAlignment,
		TsvOut:           *tsvOut,
	}

	pFaVisualize(s)
}

// pFaVisualize parses an input pFASTA file and Visualizes the file according to user-defined settings.
func pFaVisualize(s VisualizeSettings) {
	if s.TsvOut {
		browser.PFaVisualizerTsv(s.InFile, s.OutDir, s.Start, s.End, s.StartOfAlignment, s.EndOfAlignment, s.SigFigs, s.DecimalPlaces, s.LineLength, s.Chrom)
	} else {
		browser.PFaVisualizer(s.InFile, s.OutDir, s.Start, s.End, s.StartOfAlignment, s.EndOfAlignment, s.SigFigs, s.DecimalPlaces, s.LineLength, s.Chrom)
	}
}
