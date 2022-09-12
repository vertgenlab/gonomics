// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/browser"
	"github.com/vertgenlab/gonomics/common"
	"log"
)

func multFaVisualizer(infile string, outfile string, start int, end int, noMask bool, lineLength int, endOfAlignment bool) {
	browser.MultiFaVisualizer(infile, outfile, start, end, noMask, lineLength, endOfAlignment)
}

func usage() {
	fmt.Print(
		"multFaVisualizer - Provides human-readable multiple alignment from a given multiFa.\n" +
			"Keyword 'END' for the end argument makes a visualization until the end of the fasta.\n" +
			"Usage:\n" +
			"multFaVisualizer mult.fa out.txt start end\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var noMask *bool = flag.Bool("noMask", false, "Converts all bases to upper case.")
	var lineLength *int = flag.Int("lineLength", 100, "Sets to length of each alignment line.")
	var endOfAlignment bool = false
	var endPos int
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
	start := common.StringToInt(flag.Arg(2))
	if flag.Arg(3) == "END" || flag.Arg(3) == "end" || flag.Arg(3) == "End" {
		endOfAlignment = true
	} else {
		endPos = common.StringToInt(flag.Arg(3))
	}

	multFaVisualizer(infile, outfile, start, endPos, *noMask, *lineLength, endOfAlignment)
}
