package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/common"
	"log"
)

func bedDistanceFromEnds(inFile string, chromFile string, outFile string) {
	records := bed.Read(inFile)
	ref := chromInfo.ReadToMap(chromFile)
	var lengthFromEnd int64

	for i := 0; i < len(records); i++ {
		lengthFromEnd = ref[records[i].Chrom].Size - records[i].ChromEnd
		records[i].Score = common.MinInt64(lengthFromEnd, records[i].ChromStart)
	}
	bed.Write(outFile, records, 5)
}

func usage() {
	fmt.Print(
		"bedDistanceFromEnds - Returns a bed file with the score containing the distance from the end of the chromosome.\n" +
			"Usage:\n" +
			"bedDistanceFromEnds input.bed reference.chrom.sizes output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	chromFile := flag.Arg(1)
	outFile := flag.Arg(2)

	bedDistanceFromEnds(inFile, chromFile, outFile)
}
