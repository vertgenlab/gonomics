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
	}

	starrSeq.Mkref(s)
}
