package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/alleles"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func usage() {
	fmt.Print(
	"callVariants - Inputs a directory of allele count files and outputs a concatenated file that can be used as input for variant calling.\n" +
		"Usage:\n" +
		" callVariants [options] inputDirectory/ \n" +
		"options:\n")
	flag.PrintDefaults()
}

func callVariants(inDirectory string, outFile string, sigThreshold float64, afThreshold float64) {
	log.Printf("# Merging Samples\n")
	SampleMap := alleles.CreateBatchSampleMap(inDirectory)
	Variants := alleles.ScoreVariants(SampleMap, sigThreshold, afThreshold)
	if outFile == "stdout" {
		vcf.PrintVcf(Variants)
	} else {
		vcf.Write(outFile, Variants)
	}
}

func main() {
	var expectedNumArgs int = 1
	var outFile *string = flag.String("out", "stdout", "Write output to a file")
	var sigThreshold *float64 = flag.Float64("p", 0.05, "Do not output variants with p value greater than this value")
	var afThreshold *float64 = flag.Float64("af", 0.01, "Do not output variants with allele frequency less than this value")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	flag.Parse()

	inDirectory := flag.Arg(0)

	callVariants(inDirectory, *outFile, *sigThreshold, *afThreshold)

}
