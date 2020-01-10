package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/alleles"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)


func countAlleles(inFile string, outFile string, refFile string, coverageThreshold int, minMapQ int) {

	var datavcf []*vcf.Vcf

	if refFile == "" {
		usage()
		fmt.Printf("ERROR: No reference provided\n")
	} else {
		data := alleles.CountAlleles(refFile, inFile, int64(minMapQ))
		alleles.FilterAlleles(data, int32(coverageThreshold))
		datavcf = alleles.AllelesToVcf(data)
	}

	if outFile == "stdout" {
		vcf.PrintVcf(datavcf)
	} else {vcf.Write(outFile, datavcf)}
}


func usage() {
	fmt.Print(
		"countAlleles - Returns a vcf of positions with allele counts filtered by coverage and mapping quality.\n" +
			"Usage:\n" +
			" countAlleles [options] input.sam \n" +
			"options:\n")
	flag.PrintDefaults()
}


func main() {
	var expectedNumArgs int=1
	var ref *string = flag.String("ref", "", "Reference sequence in fasta format.")
	var outFile *string = flag.String("out", "stdout", "Write output to a file.")
	var minCoverage *int = flag.Int("cov", 1, "Only report positions with coverage greater than this value.")
	var minMapQ *int = flag.Int("mapQ", 20, "Only include reads with mapping quality greater than this value.")
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