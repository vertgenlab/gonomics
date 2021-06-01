// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/browser"
	"github.com/vertgenlab/gonomics/common"
	"log"
)

func multFaVisualizer(infile string, outfile string, start int, end int, noMask bool, lineLength int) {
	browser.MultiFaVisualizer(infile, outfile, start, end, noMask, lineLength)
}

//TODO: reorder this command line to be inputs first and outputs second
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
	var lineLength *int = flag.Int("lineLength", 100, "Sets to length of each alignment line.")
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
	end := common.StringToInt(flag.Arg(3))

	multFaVisualizer(infile, outfile, start, end, *noMask, *lineLength)
}
