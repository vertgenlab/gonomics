package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/alleles"
	"log"
)


func usage() {
	fmt.Print(
		"batchAlleles - Inputs a directory of allele count files and outputs a concatenated file that can be used as input for variant calling.\n" +
			"Usage:\n" +
			" batchAlleles [options] inputDirectory/ \n" +
			"options:\n")
	flag.PrintDefaults()
}

func batchAlleles(inDirectory string, outFile string, sigThreshold float64, callVar bool) {
	fmt.Printf("# Merging Samples\n")
	SampleMap := alleles.CreateSampleMap(inDirectory)

	if callVar == false {
		alleles.WriteBatchAlleleCounts(SampleMap, outFile)
	}

	if callVar == true {
		Variants := alleles.ScoreVariants(SampleMap, sigThreshold)
		alleles.WriteVarMap(Variants, outFile)
	}
}


func main() {
	var expectedNumArgs int=1
	var outFile *string = flag.String("out", "stdout", "Write output to a file")
	var sigThreshold *float64 = flag.Float64("p", 0.05, "Do not output variants with p value greater than this value")
	var callVar *bool = flag.Bool("callvars", true, "If true, p values are generated for variants using fishers exact test")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	flag.Parse()

	inDirectory := flag.Arg(0)


	batchAlleles(inDirectory, *outFile, *sigThreshold, *callVar)

}
