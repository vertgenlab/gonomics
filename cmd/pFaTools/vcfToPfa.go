package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"log"
	"os"
)

// VcfToPfaSettings defines the usage settings for the pFa VcfToPfa subcommand.
type VcfToPfaSettings struct {
	InFile    string
	RefFile	  string
	OutDir    string
	Start     int
	End       int
}

// VcfToPfaUsage defines the usage statement for the pFaTools VcfToPfa subcommand.
func VcfToPfaUsage(VcfToPfaFlags *flag.FlagSet) {
	fmt.Printf("pFaTools VcfToPfa - Builds a one-hot pFasta representing an input VCF and reference Fasta.\n" +
		"Usage:\n" +
		"PFaTools vcfToPfa inFile.vcf ref.fa outDir.pfa\n" +
		"options:\n")
	VcfToPfaFlags.PrintDefaults()
}

// parseVcfToPfaArgs is the main function of the pFaTools VcfToPfa subcommand. It parses options and runs the pFaVcfToPfa function.
func parseVcfToPfaArgs() {
	var expectedNumArgs int = 3
	var err error
	VcfToPfaFlags := flag.NewFlagSet("VcfToPfa", flag.ExitOnError)
	var start *int = VcfToPfaFlags.Int("start", 0, "Specify the position of the input sequence to begin build the pfa from. Defaults to 0.")
	var end *int = VcfToPfaFlags.Int("end", -1, "Specify the position of the input sequence to end building the pfa at. Defaults to the end.")

	err = VcfToPfaFlags.Parse(os.Args[4:])
	exception.PanicOnErr(err)
	VcfToPfaFlags.Usage = func() { VcfToPfaUsage(VcfToPfaFlags) }

	if len(VcfToPfaFlags.Args()) != expectedNumArgs {
		VcfToPfaFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(VcfToPfaFlags.Args()))
	}

	inFile := VcfToPfaFlags.Arg(0)
	refFile := VcfToPfaFlags.Arg(1)
	outDir := VcfToPfaFlags.Arg(2)

	s := VcfToPfaSettings{
		InFile:    inFile,
		RefFile:   refFile,
		OutDir:    outDir,
		Start:     *start,
		End:       *end,
	}

	vcfToPfa(s)
}

// pFaVcfToPfa parses an input pFASTA file and converts the file according to user-defined settings.
func vcfToPfa(s VcfToPfaSettings) {
	records := []pFasta.PFasta{pFasta.VcfToPfa(s.InFile, s.RefFile, s.Start, s.End)}
	pFasta.Write(s.OutDir, records)
}
