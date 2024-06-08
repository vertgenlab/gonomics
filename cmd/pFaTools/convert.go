package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"os"
	"strings"
)

// ConvertSettings defines the usage settings for the pFa Convert subcommand.
type ConvertSettings struct {
	InFile string
	OutDir string
	Start  int
	End    int
	Chrom  string
}

// ConvertUsage defines the usage statement for the pFaTools Convert subcommand.
func ConvertUsage(ConvertFlags *flag.FlagSet) {
	fmt.Printf("pFaTools Convert - Provides human-readable sequence from a given pFa.\n" +
		"Keyword 'START' for the start argument converts from the start of the pfasta.\n" +
		"Keyword 'END' for the end argument converts until the end of the pfasta.\n" +
		"Usage:\n" +
		"PFaTools convert in.fa start end outDir.txt\n" +
		"options:\n")
	ConvertFlags.PrintDefaults()
}

// parseConvertArgs is the main function of the pFaTools Convert subcommand. It parses options and runs the pFaConvert function.
func parseConvertArgs() {
	var expectedNumArgs int = 4
	var err error
	ConvertFlags := flag.NewFlagSet("Convert", flag.ExitOnError)
	var chrom *string = ConvertFlags.String("chrom", "", "Specify the name of the sequence to convert. Can be empty if only one sequence in input fasta.")
	var start int
	var end int

	err = ConvertFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	ConvertFlags.Usage = func() { ConvertUsage(ConvertFlags) }

	if len(ConvertFlags.Args()) != expectedNumArgs {
		ConvertFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(ConvertFlags.Args()))
	}

	inFile := ConvertFlags.Arg(0)
	if strings.ToLower(ConvertFlags.Arg(1)) == "start" {
		start = 0
	} else {
		start = parse.StringToInt(ConvertFlags.Arg(1))
	}
	if strings.ToLower(ConvertFlags.Arg(2)) == "end" {
		end = -1
	} else {
		end = parse.StringToInt(ConvertFlags.Arg(2))
	}
	outDir := ConvertFlags.Arg(3)

	s := ConvertSettings{
		InFile: inFile,
		OutDir: outDir,
		Start:  start,
		End:    end,
		Chrom:  *chrom,
	}

	pFaConvert(s)
}

// pFaConvert parses an input pFASTA file and converts the file according to user-defined settings.
func pFaConvert(s ConvertSettings) {
	records := []pFasta.PFasta{pFasta.faToPfa(s.InFile, s.Start, s.End, s.Chrom)}
}
