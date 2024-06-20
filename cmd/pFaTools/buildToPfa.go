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

// BuildPfaSettings defines the usage settings for the pFa BuildPfa subcommand.
type BuildPfaSettings struct {
	InFile    string
	OutDir    string
	InputType string
	Start     int
	End       int
	Chrom     string
}

// BuildPfaUsage defines the usage statement for the pFaTools BuildPfa subcommand.
func BuildPfaUsage(BuildPfaFlags *flag.FlagSet) {
	fmt.Printf("pFaTools BuildPfa - Provides human-readable sequence from a given pFa.\n" +
		"Keyword FASTA for inputType argument specifies that inFile is a Fasta.\n" +
		"Usage:\n" +
		"PFaTools buildPfa inFile outDir.pfa inputType\n" +
		"options:\n")
	BuildPfaFlags.PrintDefaults()
}

// parseBuildPfaArgs is the main function of the pFaTools BuildPfa subcommand. It parses options and runs the pFaBuildPfa function.
func parseBuildPfaArgs() {
	var expectedNumArgs int = 3
	var err error
	BuildPfaFlags := flag.NewFlagSet("BuildPfa", flag.ExitOnError)
	var start *int = BuildPfaFlags.Int("start", 0, "Specify the position of the input sequence to begin build the pfa from. Defaults to 0.")
	var end *int = BuildPfaFlags.Int("end", -1, "Specify the position of the input sequence to end building the pfa at. Defaults to the end.")
	var chrom *string = BuildPfaFlags.String("chrom", "", "Specify the name of the sequence to build. Can be empty if only one sequence in input fasta.")

	err = BuildPfaFlags.Parse(os.Args[3:])
	exception.PanicOnErr(err)
	BuildPfaFlags.Usage = func() { BuildPfaUsage(BuildPfaFlags) }

	if len(BuildPfaFlags.Args()) != expectedNumArgs {
		BuildPfaFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(BuildPfaFlags.Args()))
	}

	inFile := BuildPfaFlags.Arg(0)
	outDir := BuildPfaFlags.Arg(1)
	inputType := BuildPfaFlags.Arg(2)

	s := BuildPfaSettings{
		InFile:    inFile,
		OutDir:    outDir,
		InputType: inputType,
		Start:     *start,
		End:       *end,
		Chrom:     *chrom,
	}

	pFaBuildPfa(s)
}

// pFaBuildPfa parses an input pFASTA file and converts the file according to user-defined settings.
func pFaBuildPfa(s BuildPfaSettings) {
	if strings.ToLower(s.InputType) == "fasta" {
		records := []pFasta.PFasta{pFasta.MultiFaToPfa(s.InFile, s.Start, s.End, s.Chrom)}
		pFasta.Write(s.OutDir, records)
	} else {
		log.Fatalf("Requested format unavailable.")
	}
}
