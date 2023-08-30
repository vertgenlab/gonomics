// Command Group: "Sequence Evolution & Reconstruction"

// Returns number of mutations that separate two sequences for a given window size
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"math"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

func faFindFast(inFile string, outFile string, referenceName string, queryName string, windowSize int, refChromName string, removeN bool, longOutput bool, divergenceRate float64) {
	records := fasta.Read(inFile)

	var reference, query []dna.Base
	referenceCount := 0
	queryCount := 0

	//if reference and query names were not specified in the command, then use 1st and 2nd fasta records' names, respectively
	//will still proceed to check for non=unique names
	if referenceName == "" {
		referenceName = records[0].Name
	}
	if queryName == "" {
		queryName = records[1].Name
	}

	for i := 0; i < len(records); i++ { //for each fasta record, check if name matches reference or query
		if records[i].Name == referenceName { //if name matches reference, extract reference sequence
			if referenceCount == 1 { //however, if name has already appeared once before, fatal error for non-unique names
				log.Fatalf("Reference sequence name is not unique in the input file")
			}
			reference = records[i].Seq
			referenceCount += 1
		} else if records[i].Name == queryName {
			if queryCount == 1 {
				log.Fatalf("Query sequence name is not unique in the input file")
			}
			query = records[i].Seq
			queryCount += 1
		}
	}

	if len(reference) != len(query) {
		log.Fatalf("Reference and query sequences are not of equal length")
	}

	file := fileio.EasyCreate(outFile)
	speedyWindowDifference(windowSize, reference, query, refChromName, removeN, longOutput, divergenceRate, file)
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
	var referenceName *string = flag.String("referenceName", "", "Specify the name of the reference sequence")
	var queryName *string = flag.String("queryName", "", "Specify the name of the query sequence")
	var windowSize *int = flag.Int("windowSize", 1000, "Specify the window size")
	var refChromName *string = flag.String("chrom", "", "Specify a chrom name of the reference sequence")
	var removeN *bool = flag.Bool("removeN", false, "Excludes bed regions with Ns in the reference from the output.")
	var longOutput *bool = flag.Bool("longOutput", false, "Print percent diverged and raw -Log10PValue in output. Requires the 'divergenceRate' argument.")
	var divergenceRate *float64 = flag.Float64("divergenceRate", math.MaxFloat64, "Set the null divergence rate for p value calculations with 'longOutput'.")

	if *longOutput && *divergenceRate == math.MaxFloat64 {
		log.Fatalf("Error: must set a 'divergenceRate' if using the 'longOutput' option.")
	}

	if *divergenceRate != math.MaxFloat64 {
		if *divergenceRate < 0 || *divergenceRate > 1 {
			log.Fatalf("Error: divergence rate must be a value between 0 and 1.")
		}
	}

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

	faFindFast(inFile, outFile, *referenceName, *queryName, *windowSize, *refChromName, *removeN, *longOutput, *divergenceRate)
}
