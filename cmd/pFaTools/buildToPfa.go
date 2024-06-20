package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"log"
	"os"
	"strings"
)

// BuildToPfaSettings defines the usage settings for the pFa BuildToPfa subcommand.
type BuildToPfaSettings struct {
	InFile    string
	OutDir    string
	InputType string
	Start     int
	End       int
	Chrom     string
}

// BuildToPfaUsage defines the usage statement for the pFaTools BuildToPfa subcommand.
func BuildToPfaUsage(BuildToPfaFlags *flag.FlagSet) {
	fmt.Printf("pFaTools BuildToPfa - Provides human-readable sequence from a given pFa.\n" +
		"Keyword FASTA for inputType argument specifies that inFile is a Fasta.\n" +
		"Usage:\n" +
		"PFaTools buildToPfa inFile outDir.pfa inputType\n" +
		"options:\n")
	BuildToPfaFlags.PrintDefaults()
}

// parseBuildToPfaArgs is the main function of the pFaTools BuildToPfa subcommand. It parses options and runs the pFaBuildToPfa function.
func parseBuildToPfaArgs() {
	var expectedNumArgs int = 3
	var err error
	BuildToPfaFlags := flag.NewFlagSet("BuildToPfa", flag.ExitOnError)
	var start *int = BuildToPfaFlags.Int("start", 0, "Specify the position of the input sequence to begin build the pfa from. Defaults to 0.")
	var end *int = BuildToPfaFlags.Int("end", -1, "Specify the position of the input sequence to end building the pfa at. Defaults to the end.")
	var chrom *string = BuildToPfaFlags.String("chrom", "", "Specify the name of the sequence to build. Can be empty if only one sequence in input fasta.")

	err = BuildToPfaFlags.Parse(os.Args[3:])
	exception.PanicOnErr(err)
	BuildToPfaFlags.Usage = func() { BuildToPfaUsage(BuildToPfaFlags) }

	if len(BuildToPfaFlags.Args()) != expectedNumArgs {
		BuildToPfaFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(BuildToPfaFlags.Args()))
	}

	inFile := BuildToPfaFlags.Arg(0)
	outDir := BuildToPfaFlags.Arg(1)
	inputType := BuildToPfaFlags.Arg(2)

	s := BuildToPfaSettings{
		InFile:    inFile,
		OutDir:    outDir,
		InputType: inputType,
		Start:     *start,
		End:       *end,
		Chrom:     *chrom,
	}

	pFaBuildToPfa(s)
}

// pFaBuildToPfa parses an input pFASTA file and converts the file according to user-defined settings.
func pFaBuildToPfa(s BuildToPfaSettings) {
	if strings.ToLower(s.InputType) == "fasta" {
		records := []pFasta.PFasta{pFasta.MultiFaToPfa(s.InFile, s.Start, s.End, s.Chrom)}
		pFasta.Write(s.OutDir, records)
	} else {
		log.Fatalf("Requested format unavailable.")
	}
}
