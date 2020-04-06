package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/browser"
	"log"
)

func multFaVisualizeBeds(bedFile string, alnFile string, outFormat bool, noMask bool, lineLength int64) {
	b := bed.Read(bedFile)
	var outFile string

	for i := 0; i < len(b); i++ {
		if outFormat {
			outFile = b[i].Name + ".txt"
		} else {
			outFile = fmt.Sprintf("%s_%d_%d.txt", b[i].Chrom, b[i].ChromStart, b[i].ChromEnd)
		}
		browser.MultiFaVisualizer(alnFile, outFile, b[i].ChromStart, b[i].ChromEnd, noMask, lineLength)
	}
}

func usage() {
	fmt.Print(
		"multFaVisualizer - Provides human-readable multiple alignments for all entries in a bed file.\n" +
			"All bed entries must be on the same chromosome to interface with multiFa file.\n" +
			"Usage:\n" +
			"multiFaVisualizeBeds in.bed mult.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var noMask *bool = flag.Bool("noMask", false, "Converts all bases to upper case.")
	var outFormat *bool = flag.Bool("outFormatName", false, "Uses the name column as the outfile name (name.txt).")
	var lineLength *int64 = flag.Int64("lineLength", 100, "Sets to length of each alignment line.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	bedFile := flag.Arg(0)
	alnFile := flag.Arg(1)

	multFaVisualizeBeds(bedFile, alnFile, *outFormat, *noMask, *lineLength)
}
