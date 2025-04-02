package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/starrSeq"
	"log"
	"os"
)

func inputSeqUsage(inputSeqFlags *flag.FlagSet) {
	fmt.Print("inputSeq -- a tool to analyze input library sequencing data\n" +
		"User specifies an input sam/bam alignment file and a list of bed regions corresponsing to" +
		"constructs to be counted. An output text file that has costruct name, reads per 500bp, percentage " +
		"in library and normalization factor will be created\n" +
		"Usage:\n" +
		"inputSeq [options] in.sam in.bed out.file\n\n")
	inputSeqFlags.PrintDefaults()
}

func parseInputSeqArgs() {
	var expectedNumArgs int = 3
	inputSeqFlags := flag.NewFlagSet("inputSeq", flag.ExitOnError)
	var plasmidUMI *string = inputSeqFlags.String("plasmidUMI", "", "Provide a file name where plasmid UMI statistics will be written. "+
		"Currently, the program searches for the UMI after this sequence: GCATGCGGAT")
	var dualBx *bool = inputSeqFlags.Bool("dualBx", false, "the bed file provided for counts in a bed file of barcodes, with barcodes on each side of the construct.")

	err := inputSeqFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)

	inputSeqFlags.Usage = func() { inputSeqUsage(inputSeqFlags) }

	if len(inputSeqFlags.Args()) != expectedNumArgs {
		inputSeqFlags.Usage()
		log.Fatalf("Error: Expected %d arguments but got %d\n", expectedNumArgs, len(inputSeqFlags.Args()))
	}

	if *plasmidUMI != "" {
		starrSeq.PlasmidUMI(inputSeqFlags.Arg(0), inputSeqFlags.Arg(1), *plasmidUMI)
	}

	s := starrSeq.InputSeqSettings{
		InSam:   inputSeqFlags.Arg(0),
		InBed:   inputSeqFlags.Arg(1),
		Outfile: inputSeqFlags.Arg(2),
		DualBx:  *dualBx,
	}
	starrSeq.ParseInputSequencingSam(s)
}
