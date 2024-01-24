package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"os"
)

// ExtractSettings defines the usage settings for the pFa extract subcommand.
type ExtractSettings struct {
	InFile  string
	Chrom   string
	OutName string
	OutFile string
	Start   int
	End     int
}

// extractUsage defines the usage statement for the pFaTools extract subcommand.
func extractUsage(extractFlags *flag.FlagSet) {
	fmt.Printf("pFaTools extract - a tool for extracting a subsequence from a pFa.\n" +
		"A note on input and output file types:\n" +
		"\tpFa (probabilistic FASTA) encodes position-wise probabilities in matrix entries. Columns should therefore sum to 1.\n" +
		"Usage:\n" +
		"pFaTools extract in.pFa chrom out.fasta start end\n" +
		"options:\n")
	extractFlags.PrintDefaults()
}

// parseExtractArgs is the main function of the pFaTools extract subcommand. It parses options and runs the pFaExtract function.
func parseExtractArgs() {
	var expectedNumArgs int = 5
	var err error
	extractFlags := flag.NewFlagSet("extract", flag.ExitOnError)
	var outName *string = extractFlags.String("outName", "", "Specify name of the out file. If none is provided, output name will be the same as the input name.")
	err = extractFlags.Parse(os.Args[5:])
	exception.PanicOnErr(err)
	extractFlags.Usage = func() { extractUsage(extractFlags) }

	if len(extractFlags.Args()) != expectedNumArgs {
		extractFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(extractFlags.Args()))
	}

	inFile := extractFlags.Arg(0)
	chrom := extractFlags.Arg(1)
	start := extractFlags.Arg(2)
	end := extractFlags.Arg(3)
	outFile := extractFlags.Arg(4)

	s := ExtractSettings{
		InFile:  inFile,
		Chrom:   chrom,
		OutName: *outName,
		OutFile: outFile,
		Start:   parse.StringToInt(start),
		End:     parse.StringToInt(end),
	}

	pFaExtract(s)
}

// pFaExtract parses an input pFasta file and extracts the file according to user-defined settings.
func pFaExtract(s ExtractSettings) {
	records := []pFasta.PFasta{pFasta.Extract(pFasta.Read(s.InFile)[0], s.Start, s.End, s.OutName)}
	pFasta.Write(s.OutFile, records)
}
