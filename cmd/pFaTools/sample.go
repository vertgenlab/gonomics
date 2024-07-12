package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
)

// SampleSettings defines the usage settings for the pFa sample subcommand.
type SampleSettings struct {
	InFile     string
	Chrom      string
	OutDir     string
	NumSamples int
	SetSeed    int64
}

// sampleUsage defines the usage statement for the pFaTools sample subcommand.
func sampleUsage(sampleFlags *flag.FlagSet) {
	fmt.Printf("pFaTools sample - a tool for sampling a Fasta from a pFasta.\n" +
		"A note on input and output file types:\n" +
		"\tpFa (probabilistic FASTA) encodes position-wise probabilities in matrix entries. Columns should therefore sum to 1.\n" +
		"\tFa (FASTA) encodes integer counts in matrix entries, typically corresponding to read counts from sequencing.\n" +
		"Usage:\n" +
		"pFaTools sample in.pFa chrom out.fasta\n" +
		"options:\n")
	sampleFlags.PrintDefaults()
}

// parseSampleArgs is the main function of the pFaTools sample subcommand. It parses options and runs the pFaSample function.
func parseSampleArgs() {
	var expectedNumArgs int = 3
	var err error
	sampleFlags := flag.NewFlagSet("sample", flag.ExitOnError)
	var numSamples *int = sampleFlags.Int("numSamples", 1, "Specify the number of samples desired.")
	var setSeed *int64 = sampleFlags.Int64("setseed", 0, "Specify the seed of the random number generator.")
	err = sampleFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	sampleFlags.Usage = func() { sampleUsage(sampleFlags) }

	if len(sampleFlags.Args()) != expectedNumArgs {
		sampleFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(sampleFlags.Args()))
	}

	inFile := sampleFlags.Arg(0)
	chrom := sampleFlags.Arg(1)
	outDir := sampleFlags.Arg(2)

	s := SampleSettings{
		InFile:     inFile,
		Chrom:      chrom,
		OutDir:     outDir,
		NumSamples: *numSamples,
		SetSeed:    *setSeed,
	}

	pFaSample(s)
}

// pFaSample parses an input pFASTA file and samples the file according to user-defined settings.
func pFaSample(s SampleSettings) {
	var err error
	seed := rand.New(rand.NewSource(s.SetSeed))
	for currSample := 0; currSample < s.NumSamples; currSample++ {
		var outName = fmt.Sprintf("%s/sample_%v.fa", s.OutDir, currSample)
		records := pFasta.Sample(pFasta.Read(s.InFile), s.Chrom, seed)
		out := fileio.EasyCreate(outName)
		fasta.WriteFasta(out, records, 50)

		err = out.Close()
		exception.PanicOnErr(err)
	}
}
