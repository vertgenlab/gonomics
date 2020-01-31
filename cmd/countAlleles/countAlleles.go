package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/alleles"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func countAlleles(inFile string, outFile string, refFile string, coverageThreshold int, minMapQ int, concurrentThreads int) {

	var datavcf []*vcf.Vcf
	var data alleles.SampleMap

	if refFile == "" {
		usage()
		log.Fatalf("ERROR: No reference provided\n")
	}
	log.Printf("Counting Alleles")

	if concurrentThreads == 1 {
		data = alleles.CountAlleles(refFile, inFile, int64(minMapQ))
	} else if concurrentThreads > 1 {
		data = alleles.GoCountAlleles(refFile, inFile, int64(minMapQ), concurrentThreads)
	} else {
		log.Fatalf("Error: Requires at least 1 thread")
	}

	log.Printf("Filtering Alleles")
	alleles.FilterAlleles(data, int32(coverageThreshold))
	log.Printf("Converting to VCF")
	datavcf = alleles.AllelesToVcf(data)


	if outFile == "stdout" {
		vcf.PrintVcf(datavcf)
	} else {
		log.Printf("Writing VCF")
		vcf.Write(outFile, datavcf)
	}
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
	var expectedNumArgs int = 1
	var ref *string = flag.String("ref", "", "Reference sequence in fasta format.")
	var outFile *string = flag.String("out", "stdout", "Write output to a file.")
	var minCoverage *int = flag.Int("cov", 1, "Only report positions with coverage greater than this value.")
	var minMapQ *int = flag.Int("mapQ", 20, "Only include reads with mapping quality greater than this value.")
	var concurrent *int = flag.Int("threads", 1, "number of threads used to run process. If run with >1 thread, the process will run using concurrent goroutines")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("ERROR: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)

	countAlleles(inFile, *outFile, *ref, *minCoverage, *minMapQ, *concurrent)

}
