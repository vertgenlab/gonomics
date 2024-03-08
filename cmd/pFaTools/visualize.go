package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/browser/pFaVisualizer"
	"log"
)

// VisualizeSettings defines the usage settings for the pFa Visualize subcommand.
type VisualizeSettings struct {
	InFile     string
	OutDir     string
	Start	   int
	End	       int
	SigFigs    int
	DecimalPlaces int  
	LineLength int
	Chrom	   string
	StartofAlignment bool
	EndOfAlignment bool
}

// VisualizeUsage defines the usage statement for the pFaTools Visualize subcommand.
func VisualizeUsage(VisualizeFlags *flag.FlagSet) {
	fmt.Printf("pFaTools Visualize - Provides human-readable SOMETHING HERE from a given pFa.\n" +
		"Keyword 'START' for the end argument makes a visualization beginning from the start of the pfasta.\n" +
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
	var sigFigs int = flag.Int("sigFigs", 0, "Specify the number of significant figures to round to.")
	var decimalPlaces int = flag.Int("decimal", 7, "Specify the number of decimal places to round to.")
	var lineLength int = flag.Int("lineLength", 50, "Sets to length of each alignment line.")
	var chrom string = flag.String("chrom", "", "Specify the name of the sequence to display. Required if multiple pfastas in file.")
	var startOfAlignment bool = false
	var endOfAlignment bool = false
	var startPos int
	var endPos int

	err = VisualizeFlags.Parse(os.Args[4:])
	exception.PanicOnErr(err)
	VisualizeFlags.Usage = func() { VisualizeUsage(VisualizeFlags) }

	if len(VisualizeFlags.Args()) != expectedNumArgs {
		VisualizeFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(sVisualizeFlags.Args()))
	}

	inFile := VisualizeFlags.Arg(0)
	if strings.ToLower(Visualizeflags.Arg(1)) == "start" {
		startOfAlignment = true
		startPos = 0
	} else {
		startPos = parse.StringToInt(Visualizeflags.Arg(1))
	}
	if strings.ToLower(Visualizeflags.Arg(2)) == "end" {
		endOfAlignment = true
		endPos = -1
	} else {
		endPos = parse.StringToInt(Visualizeflags.Arg(2))
	}
	outDir := VisualizeFlags.Arg(3)

	s := VisualizeSettings{
		InFile: inFile,
		OutDir: outDir,
		Start: startPos,
		End: endPos,
		SigFigs: sigFigs,
		DecimalPlaces: decimalPlaces,  
		LineLength: lineLength,
		Chrom: chrom,
		StartOfAlignment: startOfAlignment,
		EndOfAlignment: endOfAlignment,
	}

	pFaVisualize(s)
}

// pFaVisualize parses an input pFASTA file and Visualizes the file according to user-defined settings.
func pFaVisualize(s VisualizeSettings) {
	browser.pFaVisualizer(s.Infile, s.OutDir, s.Start, s.End, s.StartOfAlignment, s.EndOfAlignment, s.SigFigs, s.DecimalPlaces, s.LineLength, s.Chrom)
}
