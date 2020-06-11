package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/browser"
	"log"
	"strconv"
)

func multFaVisualizer(infile string, outfile string, start int64, end int64, noMask bool, lineLength int64) {
	browser.MultiFaVisualizer(infile, outfile, start, end, noMask, lineLength)
}

func usage() {
	fmt.Print(
		"multFaVisualizer - Provides human-readable multiple alignment from a given multiFa.\n" +
			"Usage:\n" +
			"multFaVisualizer mult.fa out.txt start end\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var noMask *bool = flag.Bool("noMask", false, "Converts all bases to upper case.")
	var lineLength *int64 = flag.Int64("lineLength", 100, "Sets to length of each alignment line.")
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

	start, _ := strconv.ParseInt(flag.Arg(2), 10, 64)
	end, _ := strconv.ParseInt(flag.Arg(3), 10, 64)

	multFaVisualizer(infile, outfile, start, end, *noMask, *lineLength)
}
