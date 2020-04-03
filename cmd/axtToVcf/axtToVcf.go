package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

// TODO: This should be more detailed
func usage() {
	fmt.Print(
		"AXT Alignment package - b\n")
	flag.PrintDefaults()
}

// TODO: add another function so that not everything is inside main
func main() {
	var expectedNumArgs int = 1
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	snps := make([]*vcf.Vcf, 0)
	axtFile := axt.Read(flag.Arg(0))
	for i := 0; i < len(axtFile); i++ {
		snps = append(snps, axt.AxtToVcf(axtFile[i])...)
	}
	vcf.Sort(snps)
	//vcf.PrintHeader()
	vcf.PrintVcf(snps)
}
