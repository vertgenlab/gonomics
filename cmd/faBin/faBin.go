package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"path/filepath"
	"strings"
)

func faBin(genome string, path string, binNum int, minSize int) {
	records := fasta.Read(genome)
	var bins map[int][]fasta.Fasta
	bins = fasta.BinGenomeNoBreaks(records, binNum, minSize)

	for i := range bins {
		var name string
		var thisContig []fasta.Fasta
		if len(bins[i]) == 1 {
			name = bins[i][0].Name
			thisContig = bins[i]
		} else { //file name = genomeName.binNum.fa
			_, assemblyFile := filepath.Split(genome)
			assembly := strings.TrimSuffix(assemblyFile, ".fa")
			name = assembly + ".bin" + fileio.IntToString(i)
			thisContig = bins[i]
		}
		namePath := path + "/" + name + ".fa"
		fasta.Write(namePath, thisContig)
	}
}

func usage() {
	fmt.Print(
		"faBin - takes a descending size ordered fasta and returns either a specified number of fastas that are " +
			"either an entire fasta entry from the input or multiple records from the input (binNum option)" +
			" or multiple fastas which have a minimum sequence size of specified length (minSize). The number of bins " +
			"specified must be at most equal to the number of entries in the input fasta. These options cannot be combined.\n" +
			"Usage:\n" +
			"faBin <options=int> inFile.fa path\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
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

	faBin(inFile, path, *binNum, *minSize)
}
