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
	InFile  string
	RefFile string
	OutDir  string
	Start   int
	End     int
}

// VcfToPfaUsage defines the usage statement for the pFaTools VcfToPfa subcommand.
func VcfToPfaUsage(VcfToPfaFlags *flag.FlagSet) {
	fmt.Printf("pFaTools VcfToPfa - Builds a pFasta representing an input VCF and reference Fasta.\n" +
		"Note that this program assumes a sorted input VCF file.\n" +
		"Positions without VCF records are represented as one-hot pDNA vectors.\n" +
		"For positions with VCF records, we represent allele frequencies (calculated from VCF samples) in pDNA vectors.\n" +
		"Currently supports variants from a single chromosome, and will fatal if a reference fasta with multiple chromosomes is provided or if a VCF with variants from multiple chromosomes is provided.\n" +
		"Usage:\n" +
		"pFaTools vcfToPfa inFile.vcf ref.fa outDir.pfa\n" +
		"options:\n")
	VcfToPfaFlags.PrintDefaults()
}

// parseVcfToPfaArgs is the main function of the pFaTools VcfToPfa subcommand. It parses options and runs the pFaVcfToPfa function.
func parseVcfToPfaArgs() {
	var expectedNumArgs int = 3
	var err error
	VcfToPfaFlags := flag.NewFlagSet("VcfToPfa", flag.ExitOnError)
	var start *int = VcfToPfaFlags.Int("start", 0, "Specify the start position of the output pFa.")
	var end *int = VcfToPfaFlags.Int("end", -1, "Specify the end position of the output pFa.")

	err = VcfToPfaFlags.Parse(os.Args[3:])
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
		InFile:  inFile,
		RefFile: refFile,
		OutDir:  outDir,
		Start:   *start,
		End:     *end,
	}

	vcfToPfa(s)
}

// pFaVcfToPfa parses an input pFASTA file and converts the file according to user-defined settings.
func vcfToPfa(s VcfToPfaSettings) {
	records := []pFasta.PFasta{pFasta.VcfToPfa(s.InFile, s.RefFile, s.Start, s.End)}
	pFasta.Write(s.OutDir, records)
}
