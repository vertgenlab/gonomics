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
			" callVariants [options] inputDirectory/ \n" +
			"options:\n")
	flag.PrintDefaults()
}

func batchAlleles(inDirectory string, outFile string) {
	fmt.Printf("Merging Samples\n")
	SampleMap := alleles.CreateSampleMap(inDirectory)
	alleles.WriteBatchAlleleCounts(SampleMap, outFile)
}

func main() {
	var expectedNumArgs int=1
	var outFile *string = flag.String("out", "stdout", "Write output to a file")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	flag.Parse()

	inDirectory := flag.Arg(0)

	batchAlleles(inDirectory, *outFile)

}
