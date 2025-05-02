package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/starrSeq"
	"log"
	"os"
)

func makeRefUsage(makeRefFlags *flag.FlagSet) {
	fmt.Printf("starrSeqAnalysis makeRef -- make a custom starrSeq refrence assembly from constructs in library.\n" +
		"A fasta file with constructs and homology arms concatenated together, along with a bed file showing the locations of constructs will be created.\n" +
		"Usage:\n" +
		"starrSeqAnalysis makeRef upstreamHomology.fa constructs.fa downstreamHomology.fa outfilePrefix\n")
	makeRefFlags.PrintDefaults()
}

func parseMakeRefArgs() {
	var expectedNumArgs int = 4
	makeRefFlags := flag.NewFlagSet("makeRef", flag.ExitOnError)
	var sepByChrom *bool = makeRefFlags.Bool("sepByChrom", false, "Put each construct on a separate chromosome. The chromomome name will correspond to the construct")
	var dualBx *int = makeRefFlags.Int("dualBx", 0, "Create a second bed file corresponding to dual barcode locations. Assumes standard Lowe lab cloning primers are used. (20bp Fprimer, 18bp RPrimer). The input value is the length of the barcode")

	err := makeRefFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	makeRefFlags.Usage = func() { makeRefUsage(makeRefFlags) }

	if len(makeRefFlags.Args()) != expectedNumArgs {
		makeRefFlags.Usage()
		log.Fatalf("Error: Expected %d arguments but got %d\n", expectedNumArgs, len(makeRefFlags.Args()))
	}

	s := starrSeq.MakeRefSettings{
		Upstream:      makeRefFlags.Arg(0),
		Constructs:    makeRefFlags.Arg(1),
		Downstream:    makeRefFlags.Arg(2),
		OutFilePrefix: makeRefFlags.Arg(3),
		SepByChrom:    *sepByChrom,
		DualBx:        *dualBx,
	}

	starrSeq.Mkref(s)
}
