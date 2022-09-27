package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bin"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func binGenome(genome string, path string, binNum int, minSize int) {
	records := fasta.Read(genome)
	var bins map[int][]fasta.Fasta
	bins = bin.BinGenomeNoBreaks(records, binNum, minSize)

	for i := range bins {
		var name string
		var thisContig []fasta.Fasta
		if len(bins[i]) == 1 {
			name = bins[i][0].Name
			thisContig = bins[i]
		} else { //file name = first_second_...
			for j := range bins[i] {
				name = name + "_" + bins[i][j].Name
			}
			thisContig = bins[i]
		}
		namePath := path + "/" + name + ".fa"
		fasta.Write(namePath, thisContig)
	}
}

func usage() {
	fmt.Print(
		"binGenome - takes a descending size ordered fasta and returns a specified number of fastas that are " +
			"either an entire fasta entry from the input or multiple records from the input. The number of bins " +
			"specified must be at most equal to the number of entries in the input fasta. " +
			"faFilter - Returns a filtered fasta based on option parameters.\n" +
			"Usage:\n" +
			"faFilter input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var minSize *int = flag.Int("minSize", -1, "Minimum number of bases that will be in a returned fasta. Cannot be used with binNum option.")
	var binNum *int = flag.Int("binNum", 1, "Number of fasta files that will be returned containing all records of input fasta. Cannot be used with minSize option.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	path := flag.Arg(1)

	binGenome(inFile, path, *binNum, *minSize)
}
