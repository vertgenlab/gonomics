package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"log"
)

func bedFilter(infile string, outfile string, threshold *int64) {
	var records []*bed.Bed = bed.Read(infile)
	var outlist []*bed.Bed

	for i := 0; i < len(records); i++ {
		if records[i].Score >= *threshold {
			outlist = append(outlist, records[i])
		}
	}

	bed.Write(outfile, outlist, 5)
}

func usage() {
	fmt.Print(
		"bedFilter - removes bed entries below a specified threshold score.\n" +
		"Usage:\n" +
		"bedFilter input.bed output.bed\n" +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var threshold *int64 = flag.Int64("threshold", 0, "Specifies the threshold value")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)

	bedFilter(infile, outfile, threshold)
}
