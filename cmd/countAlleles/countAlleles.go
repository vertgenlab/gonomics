package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/alleles"
	"log"
)


func countAlleles(inFile string, outFile string, ref string, coverageThreshold int, minMapQ int) {

	var data alleles.SampleMap

	if ref == "" {
		usage()
		fmt.Printf("ERROR: No reference provided\n")
	} else {
		data = alleles.CountAlleles(ref, inFile, int64(minMapQ))
		alleles.FilterAlleles(data, int32(coverageThreshold))
	}


	alleles.WriteAlleleCounts(data, outFile)
}


func usage() {
	fmt.Print(
		"countAlleles - Returns a list of positions with allele counts filtered for regions with putative mutations.\n" +
			"Usage:\n" +
			" countAlleles [options] input.sam \n" +
			"options:\n")
	flag.PrintDefaults()
}


func main() {
	var expectedNumArgs int=1
	var ref *string = flag.String("ref", "", "Reference sequence in fasta format.")
	var outFile *string = flag.String("out", "stdout", "Write output to a file.")
	var minCoverage *int = flag.Int("cov", 1, "Only report positions with coverage greater than this value")
	var minMapQ *int = flag.Int("mapQ", 20, "Only include reads with mapping quality greater than this value")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("ERROR: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)

	countAlleles(inFile, *outFile, *ref, *minCoverage, *minMapQ)

}