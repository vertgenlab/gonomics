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
	fmt.Printf("")
	inputSeqFlags.PrintDefaults()
}

func parseInputSeqArgs() {
	var expectedNumArgs int = 3

	inputSeqFlags := flag.NewFlagSet("inputSeq", flag.ExitOnError)
	err := inputSeqFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	inputSeqFlags.Usage = func() { inputSeqUsage(inputSeqFlags) }

	if len(inputSeqFlags.Args()) != expectedNumArgs {
		inputSeqFlags.Usage()
		log.Fatalf("Error: Expected %d arguments but got %d\n", expectedNumArgs, len(inputSeqFlags.Args()))
	}

	s := starrSeq.InputSeqSettings{
		InSam:   inputSeqFlags.Arg(0),
		InBed:   inputSeqFlags.Arg(1),
		Outfile: inputSeqFlags.Arg(2),
	}
	starrSeq.ParseInputSequencingSam(s)
}
