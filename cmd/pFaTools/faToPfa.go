package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"log"
	"os"
)

// FaToPfaSettings defines the usage settings for the pFa FaToPfa subcommand.
type FaToPfaSettings struct {
	InFile string
	OutDir string
	Start  int
	End    int
	Chrom  string
}

// FaToPfaUsage defines the usage statement for the pFaTools FaToPfa subcommand.
func FaToPfaUsage(FaToPfaFlags *flag.FlagSet) {
	fmt.Printf("pFaTools FaToPfa - Builds a one-hot pFasta representing an input Fasta.\n" +
		"Usage:\n" +
		"PFaTools faToPfa inFile.fa outDir.pfa\n" +
		"options:\n")
	FaToPfaFlags.PrintDefaults()
}

// parseFaToPfaArgs is the main function of the pFaTools FaToPfa subcommand. It parses options and runs the pFaFaToPfa function.
func parseFaToPfaArgs() {
	var expectedNumArgs int = 3
	var err error
	FaToPfaFlags := flag.NewFlagSet("FaToPfa", flag.ExitOnError)
	var start *int = FaToPfaFlags.Int("start", 0, "Specify the position of the input sequence to begin build the pfa from. Defaults to 0.")
	var end *int = FaToPfaFlags.Int("end", -1, "Specify the position of the input sequence to end building the pfa at. Defaults to the end.")
	var chrom *string = FaToPfaFlags.String("chrom", "", "Specify the name of the sequence to build. Can be empty if only one sequence in input fasta.")

	err = FaToPfaFlags.Parse(os.Args[3:])
	exception.PanicOnErr(err)
	FaToPfaFlags.Usage = func() { FaToPfaUsage(FaToPfaFlags) }

	if len(FaToPfaFlags.Args()) != expectedNumArgs {
		FaToPfaFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(FaToPfaFlags.Args()))
	}

	inFile := FaToPfaFlags.Arg(0)
	outDir := FaToPfaFlags.Arg(1)

	s := FaToPfaSettings{
		InFile: inFile,
		OutDir: outDir,
		Start:  *start,
		End:    *end,
		Chrom:  *chrom,
	}

	faToPfa(s)
}

// faToPfa parses an input pFASTA file and converts the file according to user-defined settings.
func faToPfa(s FaToPfaSettings) {
	records := []pFasta.PFasta{pFasta.MultiFaToPfa(s.InFile, s.Start, s.End, s.Chrom)}
	pFasta.Write(s.OutDir, records)
}
