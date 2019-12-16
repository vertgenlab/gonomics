package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/alleles"
	"log"
)


func countAlleles(inFile string, outFile string, ref string) {

	var data []*alleles.AlleleCount

	if ref == "" {
		data = alleles.ExtractAlleles(inFile)
	} else {
		data = alleles.ExtractAllelesWholeGenome(inFile, ref)
	}


	if outFile == "stdout" {
		alleles.WriteVariantsToTerm(data)
	} else {
		alleles.WriteVariantsToFile(data, outFile)
	}
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
	var ref *string = flag.String("ref", "", "Reference sequence in fasta format. If given, this will generate allele counts for all regions in the genome. Can be memory intensive. If no reference is given, only positions present in the sam record are considered.")
	var outFile *string = flag.String("o", "stdout", "Write output to a file")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)

	countAlleles(inFile, *outFile, *ref)

}