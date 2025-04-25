package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/starrSeq"
	"log"
	"os"
)

func bulkOutSeqUsage(bulkOutSeqFlags *flag.FlagSet) {
	fmt.Print("")
	bulkOutSeqFlags.PrintDefaults()
}

func parseBulkOutputArgs() {
	var expectedNumArgs int = 3
	bulkOutSeqFlags := flag.NewFlagSet("bulkOutput", flag.ExitOnError)
	var normFactor *string = bulkOutSeqFlags.String("normFactor", "", "Provide the output of starrSeqAnalysis inputSeq for input normalization")
	var dualBx *bool = bulkOutSeqFlags.Bool("dualBx", false, "The bed file provided reflects barcodes that are on each side of both construct")

	err := bulkOutSeqFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)

	bulkOutSeqFlags.Usage = func() { bulkOutSeqUsage(bulkOutSeqFlags) }

	if len(bulkOutSeqFlags.Args()) != expectedNumArgs {
		bulkOutSeqFlags.Usage()
		log.Fatalf("Error: Expected %d arguments but got %d\n", expectedNumArgs, len(bulkOutSeqFlags.Args()))
	}

	s := starrSeq.BulkOutputSeqSettings{
		InSam:     bulkOutSeqFlags.Arg(0),
		InBed:     bulkOutSeqFlags.Arg(1),
		OutCounts: bulkOutSeqFlags.Arg(2),
		InputNorm: *normFactor,
		DualBx:
	}
	starrSeq.BulkOutputSeq(s)
}