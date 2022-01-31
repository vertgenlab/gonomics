// Command Group: "Sequence Evolution & Reconstruction"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func faFindFast(inFile string, outFile string, windowSize int, chromName string, removeN bool, verbose bool) {
	records := fasta.Read(inFile)

	if len(records) != 2 {
		log.Fatalf("Wrong number of sequences, expecting two, found %d.\n", len(records))
	}

	if len(records[0].Seq) != len(records[1].Seq) {
		log.Fatalf("Sequences are not of equal length")
	}

	file := fileio.EasyCreate(outFile)
	speedyWindowDifference(windowSize, records[0].Seq, records[1].Seq, chromName, removeN, verbose, file)
	err := file.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"faFindFast - Returns number of mutations that separate two sequences for a given window size\n" +
			"Usage:\n" +
			" faFindFast input.fa output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var windowSize *int = flag.Int("windowSize", 1000, "Specify the window size")
	var chromName *string = flag.String("chrom", "", "Specify the chrom name")
	var removeN *bool = flag.Bool("removeN", false, "Excludes bed regions with Ns in the reference from the output.")
	var verbose *bool = flag.Bool("verbose", false, "Enable debug prints.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	faFindFast(inFile, outFile, *windowSize, *chromName, *removeN, *verbose)
}
