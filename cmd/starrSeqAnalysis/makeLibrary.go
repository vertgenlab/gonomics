package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/starrSeq"
	"log"
	"os"
)

func makeLibraryUsage(makeLibraryFlags *flag.FlagSet) {
	fmt.Printf("starrSeqAnalysis makeLibrary [options] upstreamHomology.fa constructs.fa downstreamHomolgy.fa\n")
	makeLibraryFlags.PrintDefaults()
}

func parseMakeLibraryArgs() {
	var expectedNumArgs int = 3

	makeLibraryFlags := flag.NewFlagSet("makeLibrary", flag.ExitOnError)
	var minTm *float64 = makeLibraryFlags.Float64("minTm", -1, "Minimum acceptable melting temperature for annealing overlap")
	var maxTm *float64 = makeLibraryFlags.Float64("minTm", -1, "Maximum acceptavle melting temperature for annealing overlap")
	var captSeq *string = makeLibraryFlags.String("captSeq", "", "sequence of capture sequence to be included in oligo")
	var umiLength *int = makeLibraryFlags.Int("umiLength", 0, "Length of plasmid UMI to be included")
	var maxOligoLength *int = makeLibraryFlags.Int("maxOligoLength", 500, "Maximum oligo lenght allowed")
	err := makeLibraryFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	makeLibraryFlags.Usage = func() { makeLibraryUsage(makeLibraryFlags) }

	if len(makeLibraryFlags.Args()) != expectedNumArgs {
		makeLibraryFlags.Usage()
		log.Fatalf("Error: Expected %d arguments but got %d\n", expectedNumArgs, len(makeLibraryFlags.Args()))
	}

	s := starrSeq.MakeLibSettings{
		UpHA:        makeLibraryFlags.Arg(0),
		Constructs:  makeLibraryFlags.Arg(1),
		DownHA:      makeLibraryFlags.Arg(2),
		CS:          *captSeq,
		UmiLen:      *umiLength,
		MaxOligoLen: *maxOligoLength,
		MinTemp:     *minTm,
		MaxTemp:     *maxTm,
	}

	starrSeq.MakeLibrary(s)
}
